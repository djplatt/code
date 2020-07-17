
typedef struct{
  unsigned long int q;
  unsigned long int index;
  int_double x;
} node;

#define HEAP_LEN (1<<20)

node low_left[HEAP_LEN];
node low_right[HEAP_LEN];
node high_left[HEAP_LEN];
node high_right[HEAP_LEN];
long unsigned int llptr=0,lrptr=0,hlptr=0,hrptr=0;

void print_node(node n)
{
  printf("q: %lu, ind: %lu ",n.q,n.index);
  print_int_double_str("",n.x);
}

void print_nodel(node n)
{
  printf("%7lu %7lu %16.14e\n",n.q,n.index,n.x.left);
}

void print_noder(node n)
{
  printf("%7lu %7lu %16.14e\n",n.q,n.index,-n.x.right);
}

void print_heap(node *heap, int len)
{
  for(int i=0;i<len;i++)
    print_node(heap[i]);
}

bool llcomp(node n1,node n2)
{
  return(n1.x.left>n2.x.left);
}

bool lrcomp(node n1,node n2)
{
  return(n1.x.right<n2.x.right);
}

bool hlcomp(node n1,node n2)
{
  return(n1.x.left<n2.x.left);
}

bool hrcomp(node n1,node n2)
{
  return(n1.x.right>n2.x.right);
}


void check_lowl (node n)
{
  //printf("In check_lowl with ");print_node(n);
  if(llptr<HEAP_LEN)
    {
      low_left[llptr++]=n;
      if(llptr==HEAP_LEN)
	make_heap(low_left,&low_left[HEAP_LEN],llcomp);
    }
  else
    {
      if(llcomp(n,low_left[0]))
	{
	  //printf("Removing node      - ");print_node(low_left[0]);
	  //printf("Replacing with     - ");print_node(n);
	  pop_heap(low_left,&low_left[HEAP_LEN],llcomp);
	  low_left[HEAP_LEN-1]=n;
	  push_heap(low_left,&low_left[HEAP_LEN],llcomp);
	}
    }
}

void check_lowr (node n)
{
  if(lrptr<HEAP_LEN)
    {
      low_right[lrptr++]=n;
      if(lrptr==HEAP_LEN)
	make_heap(low_right,&low_right[HEAP_LEN],lrcomp);
    }
  else
    {
      if(lrcomp(n,low_right[0]))
	{
	  //printf("New low node found:- ");print_node(n);
	  //printf("Removing node      - ");print_node(low_right[0]);
	  pop_heap(low_right,&low_right[HEAP_LEN],lrcomp);
	  low_right[HEAP_LEN-1]=n;
	  push_heap(low_right,&low_right[HEAP_LEN],lrcomp);
	}
    }
}
void check_highl (node n)
{
  if(hlptr<HEAP_LEN)
    {
      high_left[hlptr++]=n;
      if(hlptr==HEAP_LEN)
	make_heap(high_left,&high_left[HEAP_LEN],hlcomp);
    }
  else
    {
      if(hlcomp(n,high_left[0]))
	{
	  //printf("New low node found:- ");print_node(n);
	  //printf("Removing node      - ");print_node(low_left[0]);
	  pop_heap(high_left,&high_left[HEAP_LEN],hlcomp);
	  high_left[HEAP_LEN-1]=n;
	  push_heap(high_left,&high_left[HEAP_LEN],hlcomp);
	}
    }
}

void check_highr (node n)
{
  if(hrptr<HEAP_LEN)
    {
      high_right[hrptr++]=n;
      if(hrptr==HEAP_LEN)
	make_heap(high_right,&high_right[HEAP_LEN],hrcomp);
    }
  else
    {
      if(hrcomp(n,high_right[0]))
	{
	  //printf("New low node found:- ");print_node(n);
	  //printf("Removing node      - ");print_node(low_right[0]);
	  pop_heap(high_right,&high_right[HEAP_LEN],hrcomp);
	  high_right[HEAP_LEN-1]=n;
	  push_heap(high_right,&high_right[HEAP_LEN],hrcomp);
	}
    }
}


void check_low(unsigned long int q, unsigned long int index, int_double x)
{
  node n;
  n.q=q;n.index=index;n.x=x;
  check_lowl(n);
  check_lowr(n);
}

void check_high(unsigned long int q, unsigned long int index, int_double x)
{
  node n;
  n.q=q;n.index=index;n.x=x;
  check_highl(n);
  check_highr(n);
}
