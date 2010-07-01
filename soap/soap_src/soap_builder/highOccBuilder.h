
struct highOccPatternCell {
    unsigned int l[1024*1024];
    unsigned int r[1024*1024];
    struct highOccPatternCell * next;
};


//hashtable essentials
unsigned int ho_hashTableSize=14814445;
unsigned int ho_prime=3080419483; //a relative prime to the dna length
unsigned int ho_hash_a=0, ho_hash_b=0;
unsigned int ho_ttlItem=0;
unsigned int ho_ttlOccurrence=0;

//hash func
typedef struct hashItem * hashItemPtr;
struct hashItem {
    unsigned int sa_l;
    unsigned int sa_r;
    unsigned int occCount;
    unsigned int occIndex;
    hashItemPtr next;
};


typedef struct hashTyp * hashTypPtr;
struct hashTyp {
    unsigned int count;
    unsigned int index;
    hashItemPtr item;
};

hashTypPtr ho_hashtable;


//===============================
const unsigned int ho_Size = 104857600;

typedef struct ho_cellSpace * ho_cellSpacePtr;
struct ho_cellSpace {
    unsigned int count;
    unsigned int *saL;
    unsigned int *saR;
    struct ho_cellSpace * next;
};


ho_cellSpacePtr ho_root;
ho_cellSpacePtr ho_currentRoot;

void ho_initialize() {
//    fprintf(stderr,"Initialize temporary hash.\n"); 
    srand((unsigned)time(0)); 
    ho_hash_a= rand()%ho_prime;
    ho_hash_b= rand()%ho_prime;
    
    ho_root=malloc(sizeof(struct ho_cellSpace));
    ho_root->saL = malloc(sizeof(unsigned int) * ho_Size );
    ho_root->saR = malloc(sizeof(unsigned int) * ho_Size );
    ho_root->count=0;
    ho_root->next=NULL;
    
    ho_currentRoot=ho_root;
    

    ho_ttlItem=0;
    ho_ttlOccurrence=0;
}

void ho_append(unsigned int saL, unsigned int saR) {
    if (ho_currentRoot->count >= ho_Size) {
        ho_currentRoot->next = malloc(sizeof(struct ho_cellSpace));
        ho_currentRoot=ho_currentRoot->next;
        
        ho_currentRoot->saL=malloc(sizeof(unsigned int)*ho_Size);
        ho_currentRoot->saR=malloc(sizeof(unsigned int)*ho_Size);
        ho_currentRoot->count=0;
        ho_currentRoot->next=NULL;
    }
    ho_currentRoot->saL[ho_currentRoot->count]=saL;
    ho_currentRoot->saR[ho_currentRoot->count]=saR;
    ho_currentRoot->count++;
    ho_ttlItem++;
    ho_ttlOccurrence+=saR-saL+1;
}

void ho_recurFree(ho_cellSpacePtr node) {
    if (node!=NULL) {
        ho_recurFree(node->next);
        free(node);
    }
}
void ho_recurFreeItem(hashItemPtr node) {
    if (node!=NULL) {
        ho_recurFreeItem(node->next);
        free(node);
    }
}

unsigned int ho_count() {
    return ho_ttlItem;
}
void ho_free() {
    ho_recurFree(ho_root);
    unsigned int i;
     for (i=0;i<ho_hashTableSize;i++) {
         if (ho_hashtable[i].count>0) {
            ho_recurFreeItem(ho_hashtable[i].item);
         }
     }
}



unsigned int ho_hash(unsigned int key) {
    //g(x)=(ax+b) mod p
    unsigned long long multipleA=(key * ho_hash_a)% ho_prime;
    unsigned int g=(unsigned int) (( multipleA + ho_hash_b ) % ho_prime);
    
    //f(x)=g(x) mod n
    unsigned int f=g % (ho_hashTableSize);
    
    return f;
}
