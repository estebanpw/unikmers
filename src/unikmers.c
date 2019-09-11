/*********

File        unikmers.c
Author      EPW <estebanpw@uma.es>
Description Computes HSPs using unique hits with a variable degree of uniqueness

USAGE       Usage is described by calling ./unikmers --help



**********/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <list>
#include "structs.h"
#include "alignmentFunctions.h"
#include "commonFunctions.h"
#define STARTING_SEQS 1000
#define PIECE_OF_DB_REALLOC 3200000 //half a gigabyte if divided by 8 bytes
#define RANGE 2

uint64_t custom_kmer = 32; // Defined as external in structs.h
uint64_t diffuse_z = 4; // Defined as external in structs.h


// To reduce overhead in recursion
uint64_t query_l, db_l;
char * seq_x, * seq_y;


uint64_t get_seq_len(FILE * f);
uint64_t load_seq(FILE * f, char * seq);
void save_unique_hits_pre_order(AVLTree * root, uint64_t max_unicity, FILE * out);
void save_unique_hits_in_order(AVLTree * root, uint64_t max_unicity, FILE * out);
bool compare_hits_forward(const Alignment& first, const Alignment& second);
bool compare_hits_reverse(const Alignment& first, const Alignment& second);
void transform_tree_to_hits(AVLTree * root, uint64_t max_unicity, std::list<Alignment> * list_forward, std::list<Alignment> * list_reverse);
void filter_hits(std::list<Alignment> * list_forward, std::list<Alignment> * list_reverse);
void generate_alignments_from_hits(std::list<Alignment> * list_forward, std::list<Alignment> * list_reverse, FILE * out, uint64_t min_frag_len, float min_frag_sim);
void generate_alignments_from_tree_hits(AVLTree * root, uint64_t max_unicity, FILE * out);

void reverse_complement(char * origin, char * dest);
void init_args(int argc, char ** av, FILE ** query, FILE ** database, FILE ** out_database, uint64_t * custom_kmer, uint64_t * max_unicity, uint64_t * min_frag_len, float * min_frag_sim);
void write_csv_header(FILE * out_database, uint64_t query_l, uint64_t db_l);

int main(int argc, char ** av){

    //Store positions of kmers
    uint64_t n_pools_used = 0;

    Mempool_l mp[MAX_MEM_POOLS];
    init_mem_pool_llpos(&mp[n_pools_used]);

    uint64_t n_pools_used_AVL = 0;
    Mempool_AVL mp_AVL[MAX_MEM_POOLS];
    init_mem_pool_AVL(&mp_AVL[n_pools_used_AVL]);
    
    
    AVLTree * root = NULL;


    uint64_t i, max_unicity = 2;
    float min_frag_sim = 0.75;
    uint64_t min_frag_len = 50;

    FILE * query = NULL, * database = NULL, * out_database = NULL;
    
    
    init_args(argc, av, &query, &database, &out_database, &custom_kmer, &max_unicity, &min_frag_len, &min_frag_sim);


    fprintf(stdout, "[INFO] Loading query\n");
    
    // Variables to read kmers
    char c = 'N'; //Char to read character
    
    // Read full sequence
    query_l = get_seq_len(query);

    seq_x = NULL;
    if ((seq_x = (char *) calloc(query_l, sizeof(char))) == NULL) {
        terror("Could not allocate memory for sequence x");
    }

    load_seq(query, seq_x);

    uint64_t a_hundreth = MAX(1, (query_l/100));
    unsigned char curr_kmer[custom_kmer], reverse_kmer[custom_kmer];
    curr_kmer[0] = reverse_kmer[0] = '\0';
    uint64_t word_size = 0, word_size_rev = 0;

    //To hold all information related to database
    uint64_t current_len = 0;   

    
    while( current_len < query_l ) {
        
        c = seq_x[current_len++];
            
        if(c == 'A' || c == 'C' || c == 'G' || c == 'T') {
            
            curr_kmer[word_size] = c;
            if(word_size < custom_kmer) ++word_size;
            
            if(current_len % a_hundreth == 0) { 
                fprintf(stdout, "\r%" PRIu64"%%...", 1+100*current_len/query_l); 
                fflush(stdout);
            }

        }else{ //It can be anything (including N, Y, X ...)

            if(c != '\n' && c != '>') {
                
                word_size = 0;
                ++current_len;

            } 
        }

        if(word_size == custom_kmer){
            
            //fprintf(stdout, "\t X-F %.*s -> %"PRIu64" @p:%"PRIu64"\n", 32, curr_kmer, hashOfWord(&curr_kmer[0], custom_kmer, 0), current_len);

            root = insert_AVLTree_x(root, hashOfWord(&curr_kmer[0], custom_kmer, 0), mp_AVL, &n_pools_used_AVL, current_len, mp, &n_pools_used, FORWARD);

            // Overlapping
            memmove(&curr_kmer[0], &curr_kmer[1], custom_kmer-1);
            --word_size;
        }
    }

    fprintf(stdout, "[INFO] Query computed of length %" PRIu64".\n", current_len);
    fclose(query);




    

    
    
    
    // Get file length
    
    db_l = get_seq_len(database);
    a_hundreth = MAX(1, (db_l/100));

    seq_y = NULL;
    if ((seq_y = (char *) calloc(db_l, sizeof(char))) == NULL) {
        terror("Could not allocate memory for sequence y");
    }

    load_seq(database, seq_y);

    fprintf(stdout, "[INFO] Computing database.\n");


    current_len = 0;

    curr_kmer[0] = reverse_kmer[0] = '\0';
    word_size = 0, word_size_rev = custom_kmer-1;

    
    while( current_len < db_l ) {
        
        c = seq_y[current_len++];

        
            
        if(c == 'A' || c == 'C' || c == 'G' || c == 'T') {
            
            curr_kmer[word_size] = c;
            switch(c){
                case ('A'): reverse_kmer[word_size_rev] = (unsigned)'T';
                break;
                case ('C'): reverse_kmer[word_size_rev] = (unsigned)'G';
                break;
                case ('G'): reverse_kmer[word_size_rev] = (unsigned)'C';
                break;
                case ('T'): reverse_kmer[word_size_rev] = (unsigned)'A';
                break;
            }
            if(word_size_rev != 0) --word_size_rev;

            if(word_size < custom_kmer) ++word_size;
            
            if(current_len % a_hundreth == 0) { 
                fprintf(stdout, "\r%" PRIu64"%%...", 1+100*current_len/db_l); 
                fflush(stdout);
            }

        }else{ //It can be anything (including N, Y, X ...)

            if(c != '\n' && c != '>') {
                
                word_size = 0;
                word_size_rev = custom_kmer-1;
                ++current_len;

            } 
        }

        if(word_size == custom_kmer){
            
            //fprintf(stdout, "\t Y-F %.*s -> %"PRIu64" @p:%"PRIu64"\n", 32, curr_kmer, hashOfWord(&curr_kmer[0], custom_kmer, 0), current_len);
            //fprintf(stdout, "\t Y-R %.*s -> %"PRIu64" @p:%"PRIu64"\n", 32, reverse_kmer, hashOfWord(&reverse_kmer[0], custom_kmer, 0), current_len);
            //getchar();

            root = insert_AVLTree_y(root, hashOfWord(&curr_kmer[0], custom_kmer, 0), mp_AVL, &n_pools_used_AVL, current_len, mp, &n_pools_used, FORWARD);
            root = insert_AVLTree_y(root, hashOfWord(&reverse_kmer[0], custom_kmer, 0), mp_AVL, &n_pools_used_AVL, current_len, mp, &n_pools_used, REVERSE);

            // Overlapping
            
            memmove(&curr_kmer[0], &curr_kmer[1], custom_kmer-1);
            memmove(&reverse_kmer[1], &reverse_kmer[0], custom_kmer-1);
            --word_size;
        }
    }

    fprintf(stdout, "[INFO] Database computed of length %" PRIu64".\n", current_len);
    fclose(database);

    
    write_csv_header(out_database, query_l, db_l);

    fprintf(stdout, "\n");

    /*
    fprintf(stdout, "looking for 12118931473230691009\n");
    AVLTree * wat = find_AVLTree(root, (uint64_t) 12118931473230691009L);
    if(wat != NULL){
        printf("Found it\n");
        printf("k: %"PRIu64" %.*s calc: %"PRIu64"\n", wat->key, 32, &seq_x[wat->next_in_x->pos], hashOfWord((unsigned char *) &seq_x[wat->next_in_x->pos], custom_kmer, 0));
        
    }else printf("what a pity\n");
    getchar();
    */

    //save_unique_hits_pre_order(root, max_unicity, out_database);

    std::list<Alignment> * hits_list_f = new std::list<Alignment>();
    std::list<Alignment> * hits_list_r = new std::list<Alignment>();

    transform_tree_to_hits(root, max_unicity, hits_list_f, hits_list_r);

    /*
    std::list<Alignment>::iterator it;
    std::cout << "Previous to sorted\n";
    i = 0;
    for (it=hits_list_f->begin(); it!=hits_list_f->end(); ++it) {
        std::cout << " HIT: " << (*it).x_start << " - " << (*it).y_start << "\n";
        ++i;
        if (i == 20) break;
    }
    std::cout << "\n";

    */

    // Sort
    hits_list_f->sort(compare_hits_forward);
    hits_list_r->sort(compare_hits_reverse);

    /*
    std::cout << "Post sorted\n";
    i = 1;
    for (it=hits_list_r->begin(); it!=hits_list_r->end(); ++it) {
        if (i % 100 == 0) std::cout << " HIT: " << (*it).x_start << " - " << (*it).y_start << " d: " << (int64_t) (*it).x_start + (int64_t) (*it).y_start << "\n";
        ++i;
        
    }
    std::cout << "\n";
    */
    
    fprintf(stdout, "(PRE-FILT)  Hits number: [f %" PRIu64"] [r %" PRIu64"]\n", hits_list_f->size(), hits_list_r->size());

    filter_hits(hits_list_f, hits_list_r);

    fprintf(stdout, "(POST-FILT) Hits number: [f %" PRIu64"] [r %" PRIu64"]\n", hits_list_f->size(), hits_list_r->size());

    generate_alignments_from_hits(hits_list_f, hits_list_r, out_database, min_frag_len, min_frag_sim);

    //generate_alignments_from_hits(root, max_unicity, out_database);
    fprintf(stdout, "\n");



    // Freeing

    free(seq_x);
    free(seq_y);

    for(i=0; i<n_pools_used; i++){
        free(mp[i].base);
    }
    for(i=0; i<n_pools_used_AVL; i++){
        free(mp_AVL[i].base);
    }

    




    return 0;
}

void write_csv_header(FILE * out_database, uint64_t query_l, uint64_t db_l){

    fprintf(out_database, "All by-Identity Ungapped Fragments (Hits based approach)\n");
    fprintf(out_database, "[Abr.98/Apr.2010/Dec.2011 -- <ortrelles@uma.es>\n");
    fprintf(out_database, "SeqX filename : SX\n");
    fprintf(out_database, "SeqY filename : SY\n");
    fprintf(out_database, "SeqX name : SX\n");
    fprintf(out_database, "SeqY name : SY\n");
    fprintf(out_database, "SeqX length : %" PRIu64"\n", query_l);
    fprintf(out_database, "SeqY length : %" PRIu64"\n", db_l);
    fprintf(out_database, "Min.fragment.length : unknown\n");
    fprintf(out_database, "Min.Identity : unknown\n");
    fprintf(out_database, "Tot Hits (seeds) : unknown\n");
    fprintf(out_database, "Tot Hits (seeds) used: unknown\n");
    fprintf(out_database, "Total fragments : unknown\n");
    fprintf(out_database, "========================================================\n");
    fprintf(out_database, "Total CSB: 0\n");
    fprintf(out_database, "========================================================\n");
    fprintf(out_database, "Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%%ident,SeqX,SeqY\n");
}

void save_unique_hits_pre_order(AVLTree * root, uint64_t max_unicity, FILE * out){
    unsigned char word[33];
    if(root != NULL){

        if(root->count_x != 0 && root->count_y != 0 && (root->count_x + root->count_y) <= max_unicity){
            
            perfect_hash_to_word(word, root->key, custom_kmer); word[32] = '\0';
            fprintf(stdout, "\n#H: %" PRIu64" K:%*s\n\t", root->key, (int) custom_kmer, word);
            //fprintf(stdout, "#K:%"PRIu64" ", root->key);
            llpos * aux = root->next_in_x;
            while(aux != NULL){ fprintf(stdout, "P %" PRIu64", ", aux->pos); aux = aux->next; }
            aux = root->next_in_y;
            while(aux != NULL){ fprintf(stdout, "$P %" PRIu64", ", aux->pos); aux = aux->next; }
        }
        save_unique_hits_pre_order(root->left, max_unicity, out);
        save_unique_hits_pre_order(root->right, max_unicity, out);
    }
}

void save_unique_hits_in_order(AVLTree * root, uint64_t max_unicity, FILE * out){
    unsigned char word[33];
    if(root != NULL){

        save_unique_hits_in_order(root->left, max_unicity, out);
        if(root->count_x != 0 && root->count_y != 0 && (root->count_x + root->count_y) <= max_unicity){
            perfect_hash_to_word(word, root->key, custom_kmer); word[32] = '\0';
            fprintf(stdout, "\n%*s\n", (int) custom_kmer, word);
            //fprintf(stdout, "#K:%"PRIu64" ", root->key);
            llpos * aux = root->next_in_x;
            while(aux != NULL){ fprintf(stdout, "%" PRIu64" ", aux->pos); aux = aux->next; }
            aux = root->next_in_y;
            fprintf(stdout, "\n");
            while(aux != NULL){ fprintf(stdout, "%" PRIu64" ", aux->pos); aux = aux->next; }
        }
        save_unique_hits_in_order(root->right, max_unicity, out);
    }
}


bool compare_hits_forward(const Alignment& first, const Alignment& second)
{
    int64_t d1 = (int64_t) first.x_start - (int64_t) first.y_start;
    int64_t d2 = (int64_t) second.x_start - (int64_t) second.y_start;

    if(d1 == d2) return first.x_start < second.x_start;
    return (d1 < d2);
}

bool compare_hits_reverse(const Alignment& first, const Alignment& second)
{
    int64_t d1 = (int64_t) first.x_start + (int64_t) first.y_start;
    int64_t d2 = (int64_t) second.x_start + (int64_t) second.y_start;

    if(d1 == d2) return first.x_start < second.x_start;
    return (d1 < d2);
}

void transform_tree_to_hits(AVLTree * root, uint64_t max_unicity, std::list<Alignment> * list_forward, std::list<Alignment> * list_reverse) {
    if(root != NULL){

        transform_tree_to_hits(root->left, max_unicity, list_forward, list_reverse);
        if(root->count_x != 0 && root->count_y != 0 && ( (root->count_x + root->count_y) <= max_unicity || max_unicity == 0 ) ) {
           
            llpos * aux_x = root->next_in_x;
            while(aux_x != NULL){ 
                llpos * aux_y = root->next_in_y;
                while(aux_y != NULL){ 
                    
                    Alignment a;
                    a.x_start = aux_x->pos-32;
                    a.y_start = aux_y->pos-32;
                    a.score = 0;
                    a.t_len = 0;

                    if(aux_y->strand == FORWARD) list_forward->push_back(a);
                    else list_reverse->push_back(a);

                    aux_y = aux_y->next; 
                }
                aux_x = aux_x->next; 
            }
        }
        transform_tree_to_hits(root->right, max_unicity, list_forward, list_reverse);
    }
}

void filter_hits(std::list<Alignment> * list_forward, std::list<Alignment> * list_reverse)
{

    // List forward first

    std::list<Alignment>::iterator i = list_forward->begin();
    uint64_t prev_x = (*i).x_start;
    uint64_t x;
    int64_t prev_d = (int64_t) (*i).x_start - (int64_t) (*i).y_start;
    int64_t d;
    ++i;
    while (i != list_forward->end())
    {

        d = (int64_t) (*i).x_start - (int64_t) (*i).y_start;
        x = (*i).x_start;

        if (d == prev_d && (prev_x + 32) > x)
        {
            // Remove
            list_forward->erase(i++);
        }
        else
        {
            ++i;
        }
        prev_d = d;
        prev_x = x;
    }

    // List reverse afterwards

    i = list_reverse->begin();
    prev_x = (*i).x_start;
    prev_d = (int64_t) (*i).x_start + (int64_t) (*i).y_start;
    ++i;
    while (i != list_reverse->end())
    {

        d = (int64_t) (*i).x_start + (int64_t) (*i).y_start;
        x = (*i).x_start;

        if (d == prev_d && (prev_x + 32) > x)
        {
            // Remove
            list_reverse->erase(i++);
        }
        else
        {
            ++i;
        }
        prev_d = d;
        prev_x = x;
    }

}

void generate_alignments_from_hits(std::list<Alignment> * list_forward, std::list<Alignment> * list_reverse, FILE * out, uint64_t min_frag_len, float min_frag_sim)
{

    Alignment a;

    // Forward

    std::list<Alignment>::iterator i = list_forward->begin();
    alignment_optimizer_forward((*i).x_start, (*i).y_start, query_l, db_l, seq_x, seq_y, &a);
    
    uint64_t last_reach = a.x_start + a.t_len;
    int64_t prev_d = (int64_t) (*i).x_start - (int64_t) (*i).y_start;
    int64_t d;
    
    ++i;
    while(i != list_forward->end())
    {
        d = (int64_t) (*i).x_start - (int64_t) (*i).y_start;
        
        if(d != prev_d || (d == prev_d && last_reach < (*i).x_start))
        {
            alignment_optimizer_forward((*i).x_start, (*i).y_start, query_l, db_l, seq_x, seq_y, &a);
            last_reach = a.x_start + a.t_len;
            prev_d = d;

            if(a.t_len > min_frag_len && (float)a.score / (float) a.t_len > min_frag_sim)
            {
                fprintf(out, "Frag,%" PRIu64",%" PRIu64",%" PRIu64",%" PRIu64",f,0,%" PRIu64",%" PRIu64",%" PRIu64",%f,%f,0,0\n", a.x_start, a.y_start, a.x_start+a.t_len, a.y_start+a.t_len, a.t_len, a.score, a.score, (float)a.score / (float) a.t_len, (float)a.score / (float) a.t_len);
            }
        }

        ++i;
    }

    // And reverse
    
    i = list_reverse->begin();
    alignment_optimizer_reverse((*i).x_start, (*i).y_start, query_l, db_l, seq_x, seq_y, &a);
    
    last_reach = a.x_start + a.t_len;
    prev_d = (int64_t) (*i).x_start + (int64_t) (*i).y_start;
    
    ++i;
    while(i != list_reverse->end())
    {
        d = (int64_t) (*i).x_start + (int64_t) (*i).y_start;
        
        if(d != prev_d || (d == prev_d && last_reach < (*i).x_start))
        {
            alignment_optimizer_reverse((*i).x_start, (*i).y_start, query_l, db_l, seq_x, seq_y, &a);
            last_reach = a.x_start + a.t_len;
            prev_d = d;

            if(a.t_len > min_frag_len && (float)a.score / (float) a.t_len > min_frag_sim ){
                
                fprintf(out, "Frag,%" PRIu64",%" PRIu64",%" PRIu64",%" PRIu64",r,0,%" PRIu64",%" PRIu64",%" PRIu64",%f,%f,0,0\n", a.x_start, a.y_start+a.t_len, a.x_start+a.t_len, a.y_start, a.t_len, a.score, a.score, (float)a.score / (float) a.t_len, (float)a.score / (float) a.t_len);
        
            }
        }

        ++i;
    }
}

void generate_alignments_from_tree_hits(AVLTree * root, uint64_t max_unicity, FILE * out){

    Alignment a;

    if(root != NULL){

        generate_alignments_from_tree_hits(root->left, max_unicity, out);
        if(root->count_x != 0 && root->count_y != 0 && (root->count_x + root->count_y) <= max_unicity){
           
           
            //fprintf(stdout, "On hash %"PRIu64"\n", root->key);



            llpos * aux_x = root->next_in_x;
            while(aux_x != NULL){ 
                llpos * aux_y = root->next_in_y;
                while(aux_y != NULL){ 

                    if(aux_y->strand == FORWARD) alignment_optimizer_forward(aux_x->pos-32, aux_y->pos-32, query_l, db_l, seq_x, seq_y, &a);
                    if(aux_y->strand == REVERSE) alignment_optimizer_reverse(aux_x->pos-32, aux_y->pos-32, query_l, db_l, seq_x, seq_y, &a);
                    
                    //fprintf(stdout, "ALIGNMENT: s(%"PRId64") l(%"PRIu64") x(%"PRIu64") y(%"PRIu64")\n", a.score, a.t_len, a.x_start, a.y_start);
                    //fprintf(stdout, "\t %.*s \n", (int) a.t_len, &seq_x[a.x_start]);
                    //fprintf(stdout, "\t %.*s \n", (int) a.t_len, &seq_y[a.y_start]);

                    /*
                    if(aux_y->strand == FORWARD){
                        fprintf(stdout, "@@DEBUG@@ Sweet little lies FORWARD\n");
                        fprintf(stdout, "\t %.*s -> %"PRIu64" @p:%"PRIu64"\n", 32, &seq_x[aux_x->pos - 32], hashOfWord((unsigned char *) &seq_x[aux_x->pos - 32], custom_kmer, 0), aux_x->pos);
                        fprintf(stdout, "\t %.*s -> %"PRIu64" @p:%"PRIu64"\n", 32, &seq_y[aux_y->pos - 32], hashOfWord((unsigned char *) &seq_y[aux_y->pos - 32], custom_kmer, 0), aux_y->pos);
                    }
                    */

                    /*
                    if(aux_y->strand == REVERSE){
                        fprintf(stdout, "@@DEBUG@@ Sweet little lies REVERSE\n");
                        fprintf(stdout, "\t %.*s -> %"PRIu64"\n", 32, &seq_x[aux_x->pos-32], hashOfWord((unsigned char *) &seq_x[aux_x->pos-32], custom_kmer, 0));
                        unsigned char myword[32]; myword[0] = '\0';
                        reverse_complement(&seq_y[aux_y->pos-32], (char *) myword);
                        //fprintf(stdout, "\t %.*s \n", 32, &seq_y[aux_y->pos]);
                        fprintf(stdout, "\t %.*s -> %"PRIu64"\n", 32, myword, hashOfWord(&myword[0], custom_kmer, 0));
                        getchar();
                    }
                    */
                    
                    

                    /*

                    if(aux_y->strand == REVERSE){
                        fprintf(stdout, "@@DEBUG@@ Sweet little lies REVERSE\n");
                        fprintf(stdout, "\t %.*s\n", (int) a.t_len, &seq_x[a.x_start]);
                        fprintf(stdout, "\t %.*s\n", (int) a.t_len, &seq_y[a.y_start-a.t_len]);
                        getchar();
                    }
                    */
                    
                    
                    if(a.t_len > 50 && (float)a.score / (float) a.t_len > 0.7 ){
                        if(aux_y->strand == FORWARD)
                            fprintf(out, "Frag,%" PRIu64",%" PRIu64",%" PRIu64",%" PRIu64",f,0,%" PRIu64",%" PRIu64",%" PRIu64",%f,%f,0,0\n", a.x_start, a.y_start, a.x_start+a.t_len, a.y_start+a.t_len, a.t_len, a.score, a.score, (float)a.score / (float) a.t_len, (float)a.score / (float) a.t_len);
                        else
                            fprintf(out, "Frag,%" PRIu64",%" PRIu64",%" PRIu64",%" PRIu64",r,0,%" PRIu64",%" PRIu64",%" PRIu64",%f,%f,0,0\n", a.x_start, a.y_start-a.t_len, a.x_start+a.t_len, a.y_start, a.t_len, a.score, a.score, (float)a.score / (float) a.t_len, (float)a.score / (float) a.t_len);
                    }
                    
                    aux_y = aux_y->next; 
                }
                aux_x = aux_x->next; 


            }
            
            
        }
        generate_alignments_from_tree_hits(root->right, max_unicity, out);
    }
}

void reverse_complement(char * origin, char * dest){
    int i;
    for(i=0; i<32; i++){
        if(origin[i] == 'A') dest[32-(i+1)] = 'T';
        if(origin[i] == 'C') dest[32-(i+1)] = 'G';
        if(origin[i] == 'G') dest[32-(i+1)] = 'C';
        if(origin[i] == 'T') dest[32-(i+1)] = 'A';
    }
}

void init_args(int argc, char ** av, FILE ** query, FILE ** database, FILE ** out_database, uint64_t * custom_kmer, uint64_t * max_unicity, uint64_t * min_frag_len, float * min_frag_sim){

    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           unikmers -query [query] -db [database] -out [outfile]\n");
            fprintf(stdout, "OPTIONAL:\n");
            fprintf(stdout, "           -kmer       [Integer:   k>1 (default 32)]\n");
            fprintf(stdout, "           -unique     [Integer:   u>=0 (default 2, use u=0 for no uniqueness)]\n");
            fprintf(stdout, "           -len        [Integer:   l>0 (default 50)]\n");
            fprintf(stdout, "           -sim        [Float:     0<f<1 (default 0.75)]\n");
            fprintf(stdout, "           -out        [File path]\n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            fprintf(stdout, "\n");
            fprintf(stdout, "PLEASE NOTICE: The reverse complementary is calculated for the QUERY.\n");
            exit(1);
        }
        else if(strcmp(av[pNum], "-query") == 0){
            *query = fopen64(av[pNum+1], "rt");
            if(query==NULL) terror("Could not open query file");
            pNum+=2;
        }
        else if(strcmp(av[pNum], "-db") == 0){
            *database = fopen64(av[pNum+1], "rt");
            if(database==NULL) terror("Could not open database file");
            pNum+=2;
        }
        else if(strcmp(av[pNum], "-out") == 0){
            *out_database = fopen(av[pNum+1], "wt");
            if(out_database==NULL) terror("Could not open output database file");
            pNum+=2;
        }
        else if(strcmp(av[pNum], "-kmer") == 0){
            *custom_kmer = (uint64_t) atoi(av[pNum+1]);
            if(*custom_kmer < BYTES_IN_MER) terror("K-mer size must be larger than 4");
            if(*custom_kmer % BYTES_IN_MER != 0) terror("K-mer size must be a multiple of 4");
            pNum+=2;
        }
        else if(strcmp(av[pNum], "-len") == 0){
            *min_frag_len = (uint64_t) atoi(av[pNum+1]);
            pNum+=2;
        }
        else if(strcmp(av[pNum], "-sim") == 0){
            *min_frag_sim = (uint64_t) atof(av[pNum+1]);
            pNum+=2;
        }
        else if(strcmp(av[pNum], "-unique") == 0){
            *max_unicity = (uint64_t) atoi(av[pNum+1]);
            pNum+=2;
        }
        else if(pNum > 0){
            fprintf(stderr, "Unrecognized option %s\n", av[pNum]);
            exit(-1);
        }

        if(pNum < 1) ++pNum;
        
    }
    
    if(*query==NULL || *database==NULL || *out_database==NULL) terror("A query, database and output file is required");
}

uint64_t get_seq_len(FILE * f) {
    char c = '\0';
    uint64_t l = 0;

    while(!feof(f)){
        c = getc(f);

        if(c == '>'){
            while(c != '\n') c = getc(f);
        }
        c = toupper(c);
        if(c >= 'A' && c <= 'Z'){
            ++l;
        }
    }


    rewind(f);
    return l;
}

uint64_t load_seq(FILE * f, char * seq) {
    char c = '\0';
    uint64_t l = 0;

    while(!feof(f)){
        c = getc(f);

        if(c == '>'){
            while(c != '\n') c = getc(f);
        }
        c = toupper(c);
        if(c >= 'A' && c <= 'Z'){
            seq[l] = c;
            ++l;
        }
    }


    rewind(f);
    return l;
}
