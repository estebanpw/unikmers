#define QF_LAMBDA 0.275
#define QF_KARLIN 0.333

typedef struct container{
    llpos * table[4][4][4][4][4][4][4][4][4][4][4][4];
} Container;

typedef struct index{
    Tuple_hits table[4][4][4][4][4][4][4][4][4][4][4][4];
} Index;

/*
    Nucleotides matching function
*/
int64_t compare_letters(unsigned char a, unsigned char b);

/**
 * Initialize the memory pool to later retrieve individual memory addresses for llpos
 * 
 */
void init_mem_pool_llpos(Mempool_l * mp);

/**
 * Get a new memory address from the pool mp for a type llpos
 * 
 */
llpos * getNewLocationllpos(Mempool_l * mp, uint64_t * n_pools_used);



AVLTree * getNewLocationAVLTree(Mempool_AVL * mp, uint64_t * n_pools_used, uint64_t key);


void init_mem_pool_AVL(Mempool_AVL * mp);


AVLTree * right_rotate(AVLTree * y);

AVLTree * left_rotate(AVLTree * x);

AVLTree * find_AVLTree(AVLTree * node, uint64_t key);

llpos * find_AVLTree_llpos_x(AVLTree * node, uint64_t key);

llpos * find_AVLTree_llpos_y(AVLTree * node, uint64_t key);

AVLTree * insert_AVLTree_x(AVLTree * node, uint64_t key, Mempool_AVL * mp, uint64_t * n_pools_used, uint64_t pos, Mempool_l * mp_l, uint64_t * n_pools_used_l, int strand);

AVLTree * insert_AVLTree_y(AVLTree * node, uint64_t key, Mempool_AVL * mp, uint64_t * n_pools_used, uint64_t pos, Mempool_l * mp_l, uint64_t * n_pools_used_l, int strand);

void pre_order(AVLTree * root);

uint64_t sum_of_all_tree(AVLTree * root);

void alignment_optimizer_forward(uint64_t p1, uint64_t p2, uint64_t len_x, uint64_t len_y, char * s1, char * s2, Alignment * a);

void alignment_optimizer_reverse(uint64_t p1, uint64_t p2, uint64_t len_x, uint64_t len_y, char * s1, char * s2, Alignment * a);