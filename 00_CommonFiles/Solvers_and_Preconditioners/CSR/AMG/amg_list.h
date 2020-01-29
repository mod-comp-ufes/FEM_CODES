# ifndef list_h
# define list_h

typedef struct Node {
	int elem; double val;
	struct Node *next;
} Node;

typedef struct {
	int n;
	Node *head, *tail;
} list;

typedef struct Node1 {
	int elem;
	struct Node1 *next;
} Node1;

typedef struct {
	int n;
	Node1 *head;
} list1;

typedef struct Node_p {
	double elem;
	void* val;
	struct Node_p *next;
} Node_p;

typedef struct {
	int n;
	Node_p *head, *tail;
} list_p;

list * create_l ();
void destroy_l (list *l);
void insert_l_tail (list *l, int e, double v);
void insert_l_head (list *l, int e, double v);
void insert_l_sort (list *l, int e, double v);
Node* check_in_l(list *l, int i);
list1 * create_l1 ();
void destroy_l1 (list1 *l);
void insert_l1 (list1 *l, int e);
int remove_fst (list1 *l);
list_p * create_l_p ();
void destroy_l_p (list_p* l);
void insert_l_p_tail (list_p* l, void* val, double elem);

# endif
