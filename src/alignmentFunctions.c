#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include <math.h>
#include <float.h>
#include "structs.h"
#include "alignmentFunctions.h"
#include "commonFunctions.h"



int64_t compare_letters(unsigned char a, unsigned char b){
    if(a != (unsigned char) 'N' && a != (unsigned char) '>') return (a == b) ? POINT : -POINT;
    return -POINT;
}

llpos * getNewLocationllpos(Mempool_l * mp, uint64_t * n_pools_used){

    if(mp[*n_pools_used].current == POOL_SIZE){
        *n_pools_used += 1;
        if(*n_pools_used == MAX_MEM_POOLS) terror("Reached max pools");
        init_mem_pool_llpos(&mp[*n_pools_used]);
        
    }

    llpos * new_pos = mp[*n_pools_used].base + mp[*n_pools_used].current;
    mp[*n_pools_used].current++;

    
    return new_pos;
}

void init_mem_pool_llpos(Mempool_l * mp){
    mp->base = (llpos *) calloc(POOL_SIZE, sizeof(llpos));
    if(mp->base == NULL) terror("Could not request memory pool");
    mp->current = 0;
}

AVLTree * getNewLocationAVLTree(Mempool_AVL * mp, uint64_t * n_pools_used, uint64_t key){

    if(mp[*n_pools_used].current == POOL_SIZE){
        *n_pools_used += 1;
        if(*n_pools_used == MAX_MEM_POOLS) terror("Reached max pools");
        init_mem_pool_AVL(&mp[*n_pools_used]);
        
    }

    AVLTree * new_pos = mp[*n_pools_used].base + mp[*n_pools_used].current;
    mp[*n_pools_used].current++;

    new_pos->left = NULL;
    new_pos->right = NULL;

    return new_pos;
}

void init_mem_pool_AVL(Mempool_AVL * mp){
    mp->base = (AVLTree *) calloc(POOL_SIZE, sizeof(AVLTree));
    if(mp->base == NULL) terror("Could not request memory pool");
    mp->current = 0;
}



/*
// An AVL tree node
typedef struct AVL_Node{
    uint64_t key;
    struct AVL_Node * left;
    struct AVL_Node * right;
    uint64_t height;
    llpos * next;
} AVLTree;
*/
 
// A utility function to get height of the tree

uint64_t height(AVLTree * N){
    if (N == NULL)
        return 0;
    return N->height;
}

/* Substituted by (x == NULL) ? (0) : (x->height) */
 
/* Helper function that allocates a new node with the given key and
    NULL left and right pointers. */

/* This one is substituted by AVLTree * getNewLocationAVLTree(Mempool_AVL * mp, uint64_t * n_pools_used, uint64_t key) */
 
// A utility function to right rotate subtree rooted with y
// See the diagram given above.
AVLTree * right_rotate(AVLTree * y){
    AVLTree * x = y->left;
    AVLTree * T2 = x->right;
 
    // Perform rotation
    x->right = y;
    y->left = T2;
 
    // Update heights
    //x->height = MAX((x == NULL) ? (0) : (x->left->height), (x == NULL) ? (0) : (x->right->height))+1;
    //y->height = MAX((y == NULL) ? (0) : (y->left->height), (y == NULL) ? (0) : (y->right->height))+1;
    // Update heights
    y->height = MAX(height(y->left), height(y->right))+1;
    x->height = MAX(height(x->left), height(x->right))+1;
 
    // Return new root
    return x;
}
 
// A utility function to left rotate subtree rooted with x
// See the diagram given above.
AVLTree * left_rotate(AVLTree * x){
    AVLTree * y = x->right;
    AVLTree * T2 = y->left;
 
    // Perform rotation
    y->left = x;
    x->right = T2;
 
    //  Update heights
    //x->height = MAX((x == NULL) ? (0) : (x->left->height), (x == NULL) ? (0) : (x->right->height))+1;
    //y->height = MAX((y == NULL) ? (0) : (y->left->height), (y == NULL) ? (0) : (y->right->height))+1;
    x->height = MAX(height(x->left), height(x->right))+1;
    y->height = MAX(height(y->left), height(y->right))+1;
 
    // Return new root
    return y;
}
 
// Get Balance factor of node N

int64_t get_balance(AVLTree * N){
    if (N == NULL)
        return 0;
    return height(N->left) - height(N->right);
}

/* Substituted by (node == NULL) ? (0) : ((int64_t) node->left->height - (int64_t) node->right->height) */

AVLTree * find_AVLTree(AVLTree * node, uint64_t key){
    AVLTree * found = NULL;
    if(node == NULL) return NULL;

    if (key < node->key) {
        found = find_AVLTree(node->left, key);
    } else if (key > node->key) {
        found = find_AVLTree(node->right, key);
    } else { 
        return node;
    }
    return found;
} 

llpos * find_AVLTree_llpos_x(AVLTree * node, uint64_t key){
    llpos * aux = NULL;
    if(node == NULL) return NULL;

    if (key < node->key) {
        aux = find_AVLTree_llpos_x(node->left, key);
    } else if (key > node->key) {
        aux = find_AVLTree_llpos_x(node->right, key);
    } else { 
        return node->next_in_x;
    }
    return aux;
}

llpos * find_AVLTree_llpos_y(AVLTree * node, uint64_t key){
    llpos * aux = NULL;
    if(node == NULL) return NULL;

    if (key < node->key) {
        aux = find_AVLTree_llpos_y(node->left, key);
    } else if (key > node->key) {
        aux = find_AVLTree_llpos_y(node->right, key);
    } else { 
        return node->next_in_y;
    }
    return aux;
}

// Recursive function to insert key in subtree rooted
// with node and returns new root of subtree.
AVLTree * insert_AVLTree_x(AVLTree * node, uint64_t key, Mempool_AVL * mp, uint64_t * n_pools_used, uint64_t pos, Mempool_l * mp_l, uint64_t * n_pools_used_l, int strand){
    /* 1.  Perform the normal BST insertion */
    if (node == NULL){
        
        AVLTree * n_node = getNewLocationAVLTree(mp, n_pools_used, key);
        n_node->key = key;
        n_node->count_y = 0;
        n_node->height = 1;
        n_node->count_x = 1;

        llpos * aux = getNewLocationllpos(mp_l, n_pools_used_l);
        aux->pos = pos;
        aux->strand = strand;
        n_node->next_in_x = aux;
        return n_node;
    }
 
    if (key < node->key) {
        node->left  = insert_AVLTree_x(node->left, key, mp, n_pools_used, pos, mp_l, n_pools_used_l, strand);
    } else if (key > node->key) {
        node->right = insert_AVLTree_x(node->right, key, mp, n_pools_used, pos, mp_l, n_pools_used_l, strand);
    } else { 
        
        // Equal keys are inserted as a linked list


        //printf("BINGO\n\t%"PRIu64"\n\t%"PRIu64"\n", node->key, key);

        llpos * aux = getNewLocationllpos(mp_l, n_pools_used_l);
        aux->pos = pos;
        aux->strand = strand;
        aux->next = node->next_in_x;
        node->next_in_x = aux;
        ++(node->count_x);
        
        return node;
    }
 
    

    // 2. Update height of this ancestor node 
    //node->height = 1 + MAX((node->left == NULL) ? (0) : (node->left->height), (node->right == NULL) ? (0) : (node->right->height));
    node->height = 1 + MAX(height(node->left), height(node->right));
 
    // 3. Get the balance factor of this ancestor node to check whether this node became unbalanced
    //int64_t balance = (node->left == NULL || node->right == NULL) ? (0) : ((int64_t) node->left->height - (int64_t) node->right->height);
    int64_t balance = get_balance(node);
 
    // If this node becomes unbalanced, then
    // there are 4 cases
 
    // Left Left Case
    if (balance > 1 && key < node->left->key)
        return right_rotate(node);
 
    // Right Right Case
    if (balance < -1 && key > node->right->key)
        return left_rotate(node);
 
    // Left Right Case
    if (balance > 1 && key > node->left->key)
    {
        node->left =  left_rotate(node->left);
        return right_rotate(node);
    }
 
    // Right Left Case
    if (balance < -1 && key < node->right->key)
    {
        node->right = right_rotate(node->right);
        return left_rotate(node);
    }
 
    //return the (unchanged) node pointer

    
    return node;
}


AVLTree * insert_AVLTree_y(AVLTree * node, uint64_t key, Mempool_AVL * mp, uint64_t * n_pools_used, uint64_t pos, Mempool_l * mp_l, uint64_t * n_pools_used_l, int strand){
    /* 1.  Perform the normal BST insertion */
    if (node == NULL){
        
        AVLTree * n_node = getNewLocationAVLTree(mp, n_pools_used, key);
        n_node->key = key;
        n_node->height = 1;
        n_node->count_x = 0;
        n_node->count_y = 1;

        llpos * aux = getNewLocationllpos(mp_l, n_pools_used_l);
        aux->pos = pos;
        aux->strand = strand;
        n_node->next_in_y = aux;
        return n_node;
    }
 
    if (key < node->key) {
        node->left  = insert_AVLTree_y(node->left, key, mp, n_pools_used, pos, mp_l, n_pools_used_l, strand);
    } else if (key > node->key) {
        node->right = insert_AVLTree_y(node->right, key, mp, n_pools_used, pos, mp_l, n_pools_used_l, strand);
    } else { 
        // Equal keys are inserted as a linked list
        llpos * aux = getNewLocationllpos(mp_l, n_pools_used_l);
        aux->pos = pos;
        aux->strand = strand;
        aux->next = node->next_in_y;
        node->next_in_y = aux;
        ++(node->count_y);
        
        return node;
    }
    
 
    // 2. Update height of this ancestor node
    //node->height = 1 + MAX((node->left == NULL) ? (0) : (node->left->height), (node->right == NULL) ? (0) : (node->right->height));
    node->height = 1 + MAX(height(node->left), height(node->right));
 
    // 3. Get the balance factor of this ancestor node to check whether this node became nbalanced
    //int64_t balance = (node->left == NULL || node->right == NULL) ? (0) : ((int64_t) node->left->height - (int64_t) node->right->height);
    int64_t balance = get_balance(node);
 
    // If this node becomes unbalanced, then
    // there are 4 cases
 
    // Left Left Case
    if (balance > 1 && key < node->left->key)
        return right_rotate(node);
 
    // Right Right Case
    if (balance < -1 && key > node->right->key)
        return left_rotate(node);
 
    // Left Right Case
    if (balance > 1 && key > node->left->key)
    {
        node->left =  left_rotate(node->left);
        return right_rotate(node);
    }
 
    // Right Left Case
    if (balance < -1 && key < node->right->key)
    {
        node->right = right_rotate(node->right);
        return left_rotate(node);
    }
 
    // return the (unchanged) node pointer 
    
    return node;
}
 
// A utility function to print preorder traversal
// of the tree.
// The function also prints height of every node

void pre_order(AVLTree * root){
    if(root != NULL){
        printf("#K:%" PRIu64" ", root->key);
        llpos * aux = root->next_in_x;
        while(aux != NULL){ printf("$H%" PRIu64", ", aux->pos); aux = aux->next; }
        aux = root->next_in_y;
        while(aux != NULL){ printf("$H%" PRIu64", ", aux->pos); aux = aux->next; }
        pre_order(root->left);
        pre_order(root->right);
    }
}

uint64_t sum_of_all_tree(AVLTree * root){
    uint64_t mysum = 0;
    if(root != NULL){
        
        mysum = root->count_x + root->count_y;
        mysum += sum_of_all_tree(root->left);
        mysum += sum_of_all_tree(root->right);
    }
    return mysum;
}

void alignment_optimizer_forward(uint64_t p1, uint64_t p2, uint64_t len_x, uint64_t len_y, char * s1, char * s2, Alignment * a) {

    int64_t score = 0, best_score = 0;
    int64_t curr_x = (int64_t) p1;
    int64_t curr_y = (int64_t) p2;
    int64_t best_x_start = curr_x, best_x_end = curr_x;
    int64_t best_y_start = curr_y;

    // To the right

    while(curr_x < (int64_t) len_x && curr_y < (int64_t) len_y && score >= 0 ) {
        
        if( s1[curr_x] != 'N' && s2[curr_y] != 'N' && s1[curr_x] == s2[curr_y] ) ++score; else --score;
        
        if( score >= best_score ) { best_score = score; best_x_end = curr_x; }

        ++curr_x; ++curr_y;

    }

    curr_x = (int64_t) p1;
    curr_y = (int64_t) p2;

    score = best_score;

    // To the left

    while( curr_x >= 0 && curr_y >= 0 && score >= 0 ) {
        
        if( s1[curr_x] != 'N' && s2[curr_y] != 'N' && s1[curr_x] == s2[curr_y] ) ++score; else --score;
        
        if( score >= best_score ) { best_score = score; best_x_start = curr_x; best_y_start = curr_y; }

        --curr_x; --curr_y;

    }

    a->x_start = best_x_start;
    a->y_start = best_y_start;
    a->t_len = (best_x_end - best_x_start) + 1;
    a->score = best_score;

}

void alignment_optimizer_reverse(uint64_t p1, uint64_t p2, uint64_t len_x, uint64_t len_y, char * s1, char * s2, Alignment * a) {

    int64_t score = 0, best_score = 0;
    int64_t curr_x = (int64_t) p1;
    int64_t curr_y = (int64_t) p2 + 31;
    int64_t best_x_start = curr_x, best_x_end = curr_x;
    int64_t best_y_start = curr_y;

    // To the right

    while(curr_x < (int64_t) len_x && curr_y >= 0 && score >= 0 ) {
        
        
        if( s1[curr_x] != 'N' && s2[curr_y] != 'N' && s1[curr_x] == complement(s2[curr_y]) ) ++score; else --score;
        
        if( score >= best_score ) { best_score = score; best_x_end = curr_x; best_y_start = curr_y; }

        ++curr_x; --curr_y;

    }

    curr_x = (int64_t) p1;
    curr_y = (int64_t) p2 + 31;

    // To the left
    score = best_score;

    while( curr_x >= 0 && curr_y < (int64_t) len_y && score >= 0 ) {

        if( s1[curr_x] != 'N' && s2[curr_y] != 'N' && s1[curr_x] == complement(s2[curr_y]) ) ++score; else --score;
        
        if( score >= best_score ) { best_score = score; best_x_start = curr_x; }

        --curr_x; ++curr_y;

    }

    //printf("Score achieved: %"PRId64"\n", best_score);

    a->x_start = best_x_start;
    a->y_start = best_y_start;
    a->t_len = (best_x_end - best_x_start) + 1;
    a->score = best_score;

}