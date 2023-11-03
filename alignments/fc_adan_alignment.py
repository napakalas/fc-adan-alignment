import pandas as pd
import re
import json
from tqdm import tqdm
from sentence_transformers import SentenceTransformer, util
import torch
import rdflib
from anatomy_lookup import AnatomyLookup
import ast

ADAN_FILE = 'resources/adan/adavn_vessels_with_connections.csv'
UBERON_FILE = 'resources/ontologies/uberon.ttl'
FC_ARTERY_FILE = 'resources/fc/arteries.json'

class Adan:
    def __init__(self, biobert_model):
        self.__biobert_model = biobert_model

    def __parse_vessel_name(self, vessel_name):
        # a function to split original vessel name
        # {Vessel_type}_{vessel_name}_{organ_or_tissue_region}_T{Unknown number}_{Left or right, L/R}_{vessel_number}
        # return: id, lateral, type, T, organ, name
        splitted = vessel_name.split('_')
        split_name = {}
        split_name['id'] = splitted[-1]
        split_name['lateral'] = splitted[-2]
        split_name['type'] = splitted[0]
        splitted = splitted[1:-2]
        
        # get T
        split_name['T'] = None
        if splitted[-1].startswith('T') and splitted[-1][1:].isnumeric():
            split_name['T'] = splitted[-1]
            splitted = splitted[:-1]
            
        # get organ
        pattern = r'^[A-Z]+\d+$'
        split_name['organ'] = None
        if re.match(pattern, splitted[-1]):
            organ = splitted[-1]
            splitted = splitted[:-1]
            if splitted[-1] in ['L', 'R', 'C']:
                split_name['organ'] = f'{splitted[-1]} {organ}'
                splitted = splitted[:-1]
                
        # get name
        split_name['name'] = ' '.join(splitted).lower()
        
        return list(split_name.values())

    def create_adan_index(self, adan_file=None):
        if adan_file is None: adan_file = ADAN_FILE
        df_adan = pd.read_csv(adan_file, delimiter=';')

        # split adan name and assign to dataframe
        columns = ['id', 'lateral', 'type', 'T', 'organ', 'name']
        splitted = df_adan['#[1] vessel name'].apply(lambda x: self.__parse_vessel_name(x))
        splitted = list(zip(*splitted.values))
        df_adan['id'] = splitted[0]
        df_adan['lateral'] = splitted[1]
        df_adan['type'] = splitted[2]
        df_adan['T'] = splitted[3]
        df_adan['organ'] = splitted[4]
        df_adan['name'] = splitted[5]

        # get parent and children
        def get_feature(a_list: list, extracted_col, filling_col):
            return [df_adan[df_adan[extracted_col] == element][filling_col].values[0].lower() for element in ast.literal_eval(a_list)]

        df_adan['parent name'] = df_adan['[15] input vessels'].apply(lambda names: get_feature(names, '#[1] vessel name', 'name'))
        df_adan['parent id'] = df_adan['[15] input vessels'].apply(lambda names: get_feature(names, '#[1] vessel name', 'id'))
        df_adan['children name'] = df_adan['[16] output vessels'].apply(lambda names: get_feature(names, '#[1] vessel name', 'name'))
        df_adan['children id'] = df_adan['[16] output vessels'].apply(lambda names: get_feature(names, '#[1] vessel name', 'id'))

        # flatten vessels having several children into different rows combining adan, the same row with the same 
        columns = ['#[1] vessel name', '[2] first node', '[3] second node', 'name', 'parent id', 'parent name', 'children id', 'children name', 'type']
        self.__df_adan_flatten = pd.DataFrame(columns=columns)
        for index, row in df_adan[columns].reset_index().iterrows():
            #currently, only focuss on arteries
            if row['type'] not in ['A', 'P', 'T']:
                continue
            for parent in zip(row['parent id'], row['parent name']):
                if len(row['children id']) > 0:
                    for child in zip(row['children id'], row['children name']):
                        add_row = [row['#[1] vessel name'], row['[2] first node'], row['[3] second node'], row['name'], parent[0], parent[1], child[0], child[1], row['type']]
                        self.__df_adan_flatten.loc[len(self.__df_adan_flatten.index)] = add_row
                else: # some of vessels do not have children
                    add_row = [row['#[1] vessel name'], row['[2] first node'], row['[3] second node'], row['name'], parent[0], parent[1], None, None, row['type']]
                    self.__df_adan_flatten.loc[len(self.__df_adan_flatten.index)] = add_row
                
        # create embedding of artery name in adan
        self.__adan_terms = list(set(self.__df_adan_flatten['name']).union(set(self.__df_adan_flatten['parent name'])).union(set(self.__df_adan_flatten['children name'])))
        self.__adan_terms = [term for term in self.__adan_terms if term is not None]

        #calculate embedding of adan terms
        self.__modified_adan_terms = [term if 'arter' in term else term + ' artery' for term in self.__adan_terms]
        self.__adan_term_embs = self.__biobert_model.encode(self.__modified_adan_terms, show_progress_bar=True, convert_to_tensor=True)

        # Expanding ADAN terms, by checking their availability in UBERON and then adding as additional synonyms
        g = rdflib.Graph()
        g.parse(UBERON_FILE, 'r')
        def get_synonyms(term_id):
            subject = rdflib.URIRef(term_id)
            predicate = rdflib.URIRef('http://www.geneontology.org/formats/oboInOwl#hasRelatedSynonym') | rdflib.URIRef('http://www.geneontology.org/formats/oboInOwl#hasExactSynonym')
            synonyms = []
            for o in g.objects(subject=subject, predicate=predicate):
                synonyms += [str(o)]
            return synonyms
        lookup = AnatomyLookup()
        self.__adan_synonyms_map = {}
        self.__adan_synonym_terms = []

        for term in tqdm(self.__modified_adan_terms):
            if term not in self.__adan_synonym_terms:
                self.__adan_synonyms_map[term] = term
                self.__adan_synonym_terms += [term]
                
            term_id, label, score = lookup.search(term, refine=False)
            if score > 0.97:
                synonyms = [label] + get_synonyms(term_id)
                for synonym in synonyms:
                    synonym = synonym.lower() if 'arter' in synonym.lower() else synonym.lower() + ' artery'
                    if synonym not in self.__adan_synonym_terms:
                        self.__adan_synonyms_map[synonym] = term
                        self.__adan_synonym_terms += [synonym]

        # now convert synonym to embedding
        self.__adan_synonym_term_embs = self.__biobert_model.encode(self.__adan_synonym_terms, show_progress_bar=True, convert_to_tensor=True)

        # create embeddings of with the combinations of artery's name, parent, and child
        self.__adan_triple, self.__adan_name_parent, self.__adan_parent_child, self.__adan_name_child = [], [], [], []
        self.__adan_triple_embs, self.__adan_name_parent_embs, self.__adan_parent_child_embs, self.__adan_name_child_embs = [], [], [], []
        for idx, row in self.__df_adan_flatten.iterrows():
            triple = (row['name'], row['parent name'], row['children name'])
            if triple not in self.__adan_triple:
                self.__adan_triple += [triple]
                self.__adan_triple_embs += [torch.mean(torch.stack([self.__adan_term_embs[self.__adan_terms.index(t)] for t in triple if t is not None]), 0)]
            
            name_parent = (row['name'], row['parent name'])
            if name_parent not in self.__adan_name_parent:
                self.__adan_name_parent += [name_parent]
                self.__adan_name_parent_embs += [torch.mean(torch.stack([self.__adan_term_embs[self.__adan_terms.index(t)] for t in name_parent if t is not None]), 0)]
            
            parent_child = (row['parent name'], row['children name'])
            if parent_child not in self.__adan_parent_child:
                self.__adan_parent_child += [parent_child]
                self.__adan_parent_child_embs += [torch.mean(torch.stack([self.__adan_term_embs[self.__adan_terms.index(t)] for t in parent_child if t is not None]), 0)]
                
            name_child = (row['name'], row['children name'])
            if name_child not in self.__adan_name_child:
                self.__adan_name_child += [name_child]
                self.__adan_name_child_embs += [torch.mean(torch.stack([self.__adan_term_embs[self.__adan_terms.index(t)] for t in name_child if t is not None]), 0)]
            
        self.__adan_triple_embs = torch.stack(self.__adan_triple_embs)
        self.__adan_name_parent_embs = torch.stack(self.__adan_name_parent_embs)
        self.__adan_parent_child_embs = torch.stack(self.__adan_parent_child_embs)
        self.__adan_name_child_embs = torch.stack(self.__adan_name_child_embs)

    def get_similar_adan(self, queries, terms, embs, top_k=5, min_sim=0.7):
        # a common function to get similar object to queries

        queries = [queries] if isinstance(queries, str) else queries
        queries = [query.lower() if 'arter' in query.lower() else query.lower() + ' artery' for query in queries if query is not None]
        query_embs = []
        for query in queries:
            if query is not None:
                # check in synonym
                query_emb = self.__biobert_model.encode(query, convert_to_tensor=True)
                cos_scores = util.cos_sim(query_emb, self.__adan_synonym_term_embs)[0]
                top_result = torch.topk(cos_scores, k=1)
                if top_result[0] > 0.97:
                    new_query = self.__adan_synonyms_map[self.__adan_synonym_terms[top_result[1]]]
                    query_emb = self.__adan_term_embs[self.__modified_adan_terms.index(new_query)]
                query_embs += [query_emb]
            
        query_embs = torch.mean(torch.stack(query_embs), 0)
        cos_scores = util.cos_sim(query_embs, embs)[0]
        top_results = torch.topk(cos_scores, k=top_k)

        results = []
        for score, idx in zip(top_results[0], top_results[1]):
            if score < min_sim:
                break
            results += [(terms[idx], score.item())]
        return results

    # functions to get similar terms, triples, name_parent, parent_children 
    def get_similar_adan_term(self, query, top_k=5, min_sim=0.7):
        return self.get_similar_adan(query, self.__adan_terms, self.__adan_term_embs, top_k=top_k, min_sim=min_sim)

    def get_similar_adan_triple(self, query, top_k=5, min_sim=0.7):
        return self.get_similar_adan(query, self.__adan_triple, self.__adan_triple_embs, top_k=top_k, min_sim=min_sim)

    def get_similar_adan_parent_child(self, query, top_k=5, min_sim=0.7):
        return self.get_similar_adan(self, query, self.__adan_parent_child, self.__adan_parent_child_embs, top_k=top_k, min_sim=min_sim)

    def get_similar_adan_name_parent(self, query, top_k=5, min_sim=0.7):
        return self.get_similar_adan(query, self.__adan_name_parent, self.__adan_name_parent_embs, top_k=top_k, min_sim=min_sim)

    def get_similar_adan_name_child(self, query, top_k=5, min_sim=0.7):
        return self.get_similar_adan(query, self.__adan_name_child, self.__adan_name_child_embs, top_k=top_k, min_sim=min_sim)

    def get_parent_name(self, name):
        return list(self.__df_adan_flatten[self.__df_adan_flatten['name']==name]['parent name'])[0]
