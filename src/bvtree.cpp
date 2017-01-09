/***************************************************************************
 *   DynFMI - Dynamic FM-Index for a Collection of Texts                   *
 *   Copyright (C) 2006  Wolfgang Gerlach                                  *
 *                                                                         *
 *   This program is free software: you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation, either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>. *
 ***************************************************************************/


#include "bvtree.h"


using namespace std;

//std::ostream& operator<<(std::ostream& os, node& n) { ... }

namespace bvtree {

  BVTree::~BVTree(){
	delete tempbit;
  }

	

/**
BVTree::BVTree(char *readfile){
  
  std::ifstream file(readfile);
  if (!file) 
  {
	  cout << "error reading file "<< readfile <<" !\n";
	  exit(EXIT_FAILURE);
  }


  tempnil = (BVNode*) ((RBTree*) this)->nil;
  tempbit = new bitset<2*logn>;
  

  char c;
  bitset<8> bs;

  while (file.get(c)) 
  {
	  bs = c;

	  for (int i=7; i>=0; i--) 
		  appendBit(bs[i]);
  }
  
}
**/




void BVTree::iterateReset(){
  	iterate=1;
	iterateLocal=0;
	if (root!=nil) {
		iterateNode = (BVNode*) treeMinimum(root);
		iterateRank = ((*iterateNode->block)[0])?1:0;
	} else {
		iterateNode=getNil();
		iterateRank = 0;
		}
}
  
bool BVTree::iterateGetBit(){
	return ((*iterateNode->block)[iterateLocal]);
}

size_t BVTree::iterateGetRank(bool bit){
	if (bit) return iterateRank;
	else return (iterate - iterateRank);
}

bool BVTree::iterateNext(){
	iterate++;

	if (iterate > getPositions()) return false;
	
	
	if (iterateLocal < iterateNode->myPositions-1) {
		iterateLocal++;
	} else {
		// jump to next leaf;
		iterateNode = (BVNode*) treeSuccessor(iterateNode);
		#ifndef NDEBUG
		if (iterateNode==getNil()) {
			cout << "iterateNode==getNil()" << endl;
			return false;
		}
		#endif
		iterateLocal=0;
	}
	if  ((*iterateNode->block)[iterateLocal]) iterateRank++;
	
	return true;
}

size_t BVTree::getPositions(){
	
	if (getRoot() == getNil()) return 0;
	return getRoot()->subTreePositions;
}

size_t BVTree::getRank(){
	
	if (getRoot() == getNil()) return 0;
	return getRoot()->subTreeRank;
}


void BVTree::printNode(BVNode *n){
	int commas=(2*logn)/10;
	int commashift = -1;
	
	
	cout << "address: " << n << endl;

	if (n->block != 0 && (n!=getNil()))
	{
		int size = 2*logn + 1 +commas;
		
		char *myblock = (char *) new char[size];
		if (myblock == 0) {
			cerr << "error: printNode: myblock == 0" << endl;
			exit(0);
		}
		
		//read: i, write: i+commashift
		for (int i=0; i<(int)(n->myPositions); i++){
			if (i%10 == 0) commashift++;
			myblock[i+commashift]='0' + (*n->block)[i];
			if (i+commashift >= size-1) {
				cerr << "X) printNode: array wrong index" << endl;
				exit(0);
			}
			
		}
		
		cout << "size=" << size << endl;
		cout << "n->myPositions=" << n->myPositions << endl;
		for (int i=n->myPositions; i<(2*logn); i++){
			if ((i%10) == 0) commashift++;
			
			myblock[i+commashift]='-';
			if (i+commashift > size-2) {
				cerr << "A) printNode: array wrong index" << endl;
				exit(0);
			}
		}
		
		commashift = 0;
		for (int i=10; i < 2*logn; i++){
			if (i%10 == 0) 
			{
				if (i+commashift >= size-2) {
					cerr << "B) printNode: array wrong index" << endl;
					exit(0);
				}
				
				myblock[i+commashift]=',';
				commashift++;
			}

		}

		myblock[size - 1]='\0';
		
		cout << "block: \"" << myblock << "\"" << endl;
		delete myblock;
		
	}
	else 
		cout << "block: none" << endl;
	
	
	cout << "myPositions: " << n->myPositions << endl;
	cout << "myRank: " << n->myRank << endl;
	cout << "subTreePositions: " << n->subTreePositions << endl;
	cout << "subTreeRank: " << n->subTreeRank << endl;
	cout << "color: " << ((n->color==rbtree::RED)?"RED":"BLACK") << endl;
	cout << "parent: " << n->getParent() << endl;
	cout << "left: " << n->getLeft() << endl;
	cout << "right:" << n->getRight() << endl << endl;

}


int BVTree::getTreeMaxDepth(){
	return getNodeMaxDepth(root);
}

int BVTree::getTreeMinDepth(){
	return getNodeMinDepth(root);
}


void BVTree::updateCounters(BVNode *n){

	if (n == getNil()) return;

	size_t lR = 0;
	size_t lP = 0;
	size_t rR = 0;
	size_t rP = 0;

	if (n->getLeft() != getNil()) {
		lR = n->getLeft()->subTreeRank;
		lP = n->getLeft()->subTreePositions;
	}

	if (n->getRight() != getNil()) {
		rR = n->getRight()->subTreeRank;
		rP = n->getRight()->subTreePositions;
	}

	n->subTreeRank     =lR + rR + n->myRank;
	n->subTreePositions=lP + rP + n->myPositions;
	
}

size_t BVTree::getLocalRank(BVNode* n, size_t position){
	
	#ifndef NDEBUG
	if (position > n->myPositions) {
		cerr << "error: getLocalRank: invalid position in block.\n";
		exit(EXIT_FAILURE);
	}
	#endif
	
	*tempbit =*(n->block)<<((2*logn)-position);
	return tempbit->count(); 

	// old version:
	//rank = 0;
	//for (size_t i = 0; i < position; i++) {
	//	if ((*n->block)[i]) rank++;
	//}
}

size_t BVTree::getLocalSelect1(BVNode* n, size_t query){
	size_t select =0;
	#ifndef NDEBUG
	if (query > n->myPositions) {
		cerr << "error: getLocalSelect1: invalid position in block.\n";
		exit(EXIT_FAILURE);
	}
	#endif
	size_t i;
	for (i = 0; select < query; i++) { // TODO is there a faster solution similar to rank ?
		if ((*n->block)[i]) select++;
	}

	return i;
}

size_t BVTree::getLocalSelect0(BVNode* n, size_t query){
	size_t select =0;
	#ifndef NDEBUG
	if (query > n->myPositions) {
		cerr << "error: getLocalSelect0: invalid position in block.\n";
		exit(EXIT_FAILURE);
	}
	#endif
	size_t i;
	for (i = 0; select < query; i++) {
		if (!(*n->block)[i]) select++;
	}

	return i;
}


void BVTree::printNode(size_t i){
	BVNode* x = getRoot();

	#ifndef NDEBUG	
	if (x == getNil()) { 
		cerr << "error: printNode(int i): root=NULL.\n";
		exit(EXIT_FAILURE);
		
	}
	
	if (i > x->myPositions) {
		cerr << "error: printNode(int i): invalid position in block.\n";
		exit(EXIT_FAILURE);
	}
	#endif	
	
	size_t lP=0;
	bool search = true;
	//find the corresponding block:
	while (search) 
	{
		
		lP = x->getLeft()->subTreePositions;
		

		if (lP >= i)
		{
			x=x->getLeft();
		}
		else if (lP+x->myPositions >= i){
			i-=lP;
			search = false;
		}
		else{
			i-=(lP+x->myPositions);
			x=x->getRight();
		}
	}


	cout << "i=" << i << endl;
	printNode(x);
}


bool BVTree::operator[](size_t i){

	BVNode* x = getRoot();
	size_t lsP=0;
	bool search = true;
	//find the corresponding block:
	while (search) 
	{
		
		lsP = x->getLeft()->subTreePositions;

		if (lsP >= i)
		{
			#ifndef NDEBUG
			if (x->getLeft()==getNil()) {
				cout << "lsP: " << lsP << endl;
				printNode(x);
				cout << "ihhh" << endl;
				checkTree();
				exit(EXIT_FAILURE);
			}
			#endif
			x=x->getLeft();
		}
		else if (lsP+x->myPositions >= i){
			i-=lsP;
			search = false;
		}
		else{
			i-=(lsP+x->myPositions);
			#ifndef NDEBUG
			if (x->getRight()==getNil()) {
				cout << "i: " << i << endl;
				cout << "lsP: " << lsP << endl;
				cout << "x->myPositions: " << x->myPositions << endl;
				printNode(x);
				cout << "ihhh" << endl;
				checkTree();
				exit(EXIT_FAILURE);
				}
			#endif
			x=x->getRight();
		}
	}

	return (*x->block)[i-1];  
}


size_t BVTree::rank1(size_t i){
	BVNode* x = getRoot();

	if (i == this->getPositions() + 1) i--;
	#ifndef NDEBUG	
	if (i > this->getPositions() ) {
		cerr << "error: rank1(0): invalid position in bittree: " << i << endl;
		cerr << "error: rank1(0): this->getPositions(): " <<  this->getPositions()<< endl;
		exit(EXIT_FAILURE);
	}
	#endif	


	size_t lsP=0;
	size_t lsR=0;
	size_t rank=0;
	bool search = true;
	//find the corresponding block:
	while (search) 
	{
		
		lsP = x->getLeft()->subTreePositions;
		lsR = x->getLeft()->subTreeRank;

		if (lsP >= i)
		{// cout << "L" << endl;
			x=x->getLeft();
		}
		else if (lsP+x->myPositions >= i){
			i-=lsP;
			rank+=lsR;
			search = false;
		}
		else{//cout << "R" << endl;
			i-=(lsP+x->myPositions);
			rank+=(lsR+x->myRank);
			x=x->getRight();
		}
	}

	rank+=getLocalRank(x, i);
	return rank;
}

size_t BVTree::rank0(size_t i){
	if (this->getPositions() == 0) return 0;
	return (i-rank1(i));
}

size_t BVTree::select1(size_t i){
	BVNode* x = getRoot();
	size_t select=0;
	#ifndef NDEBUG	
	if (i > x->subTreeRank ) {
		cerr << "error: select1: invalid position in bittree: " << i << endl;
		exit(EXIT_FAILURE);
	}
	#endif

	size_t lsP=0;
	size_t lsR=0;
	bool search = true;
	//find the corresponding block:
	while (search) 
	{
		// update subTree-counters
		
		
		lsP = x->getLeft()->subTreePositions;
		lsR = x->getLeft()->subTreeRank;

		if (lsR >= i) {
			x=x->getLeft();
		}
		else if (lsR+x->myRank >= i) {
			i-=lsR;
			select+=lsP;
			search = false;
		}
		else {
			i-=(lsR+x->myRank);
			select+=(lsP+x->myPositions);
			x=x->getRight();
		}
	}
	select+=getLocalSelect1(x, i);



	return select;
}

size_t BVTree::select0(size_t i){
	BVNode* x = getRoot();
	size_t select=0;

	#ifndef NDEBUG
	if (i > (x->subTreePositions - x->subTreeRank)) {
		cerr << "error: select1: invalid position in bittree: " << i << endl;
		exit(EXIT_FAILURE);
	}	
	#endif


	size_t lsP=0;
	size_t lsR=0;
	size_t lmR=0;
	size_t lmP=0;
	bool search = true;
	//find the corresponding block:
	while (search) 
	{
		// update subTree-counters
		
		
		lsP = x->getLeft()->subTreePositions;
		lsR = x->getLeft()->subTreeRank;
		lmR = x->getLeft()->myRank;
		lmP = x->getLeft()->myPositions;

		if (lsP-lsR >= i)
		{
			x=x->getLeft();
		}
		else if ((lsP-lsR)+(x->myPositions-x->myRank) >= i){
			i-=lsP-lsR;
			select+=lsP;
			search = false;
		}
		else{
			i-=((lsP-lsR)+(x->myPositions-x->myRank));
			select+=(lsP+x->myPositions);
			x=x->getRight();
		}
	}

	select+=getLocalSelect0(x, i);

	return select;	
}


void BVTree::updateCountersOnPathToRoot(BVNode *walk){

	while (walk != getNil()) {
		updateCounters(walk);
		walk=walk->getParent();
	}
}

//deletes the BVNode and all its children, destroys reb-black-tree property!
void BVTree::deleteNode(BVNode *n){
	if (n==getNil()) return;

	if (n->getLeft() != getNil()) deleteNode(n->getLeft());
	if (n->getRight() != getNil()) deleteNode(n->getRight());
	
	delete n;
}




// TODO improve by returning bitvalue

void BVTree::deleteBit(size_t i){
	size_t old_i = i; 
	BVNode* x = getRoot();
	bool bit;
	size_t rank=0;

	#ifndef NDEBUG	
	if (x == getNil()) {
		cerr << "error: deleteBit, root is empty\n"; //shouldn't happen
		exit(EXIT_FAILURE);
	}

	if (i > x->subTreePositions || i < 1) 
	{
		cerr << "error: A, position " << i <<" in block not available, only " << x->myPositions <<" positions.\n"; //shouldn't happen
		exit(EXIT_FAILURE);
	}
	#endif	

	size_t lsP=0;
	size_t lsR=0;

	bool search = true;
	//find the corresponding block:
	while (search) 
	{
		// update of pointers is not yet possible: call updateCountersOnPathToRoot

		lsP = x->getLeft()->subTreePositions;
		lsR = x->getLeft()->subTreeRank;
		

		if (lsP >= i)
		{
			#ifndef NDEBUG
			if (x->getLeft()==getNil()) exit(EXIT_FAILURE);
			#endif
			x=x->getLeft();
		}
		else if (lsP+x->myPositions >= i){
			i-=lsP;
			rank+=lsR;
			search = false;
		}
		else{
			i-=(lsP+x->myPositions);
			rank+=(lsR+x->myRank); // for speedup!
			#ifndef NDEBUG
			if (x->getRight()==getNil()) exit(EXIT_FAILURE);
			#endif
			x=x->getRight();
		}
	}


	
        // now delete the bit from the block x:
	bit =(*x->block)[i-1];

	// store bit and rank information for speedup
	lastBitDeleted=bit;
	rank+=getLocalRank(x, i);
	lastRank=(bit?rank:old_i-rank);
	
	#ifndef NDEBUG	
	if (i > x->myPositions) {
		cerr << "error: B, position " << i <<" in block not available, only " << x->myPositions <<" positions.\n"; //shouldn't happen
		exit(EXIT_FAILURE);
	}
	#endif
	bitset<2*logn> mask;
	

	if ( i > 1 ) 
	{
		mask.set();
		mask>>=(2*logn - i + 1);
		mask &= *(x->block);
		(*(x->block)>>=i)<<=i-1;  // shift bits by one	
		(*x->block)|=mask;        // restore bits in front of deleted bit
	} else {
		*(x->block)>>=1;  // shift bits by one	
	}
		
	x->myPositions--;
	if (bit) x->myRank--;

	updateCountersOnPathToRoot(x);

	
	if (x->myPositions == 0){ // if merging was not possible:
	
		rbDelete(x, callUpdateCountersOnPathToRoot);
		delete x;

	} else if (x->myPositions <= logn/2) // merge (and rotate), if necessary:
	{
		
		BVNode *sibling = (BVNode*) treeSuccessor(x);
		//cout << "try to merge -----------------------------------------" << endl;
		if (sibling != getNil())
		{
			if (x->myPositions + sibling->myPositions < 2*logn) //merge !
			{
			//cout << "merge with right sibling -----------------------------------------" << endl;
				// move bits from sibling to x:

				(*sibling->block)<<=x->myPositions;
				(*x->block)|=(*sibling->block);

				x->myPositions+=sibling->myPositions;
				x->myRank+=sibling->myRank;
				updateCountersOnPathToRoot(x);
				rbDelete(sibling, callUpdateCountersOnPathToRoot);
				delete sibling;
			}
			else sibling = getNil(); //to try the left sibling!
		}
		else if ((sibling= (BVNode*) treePredecessor(x)) != getNil())
		{
			if (x->myPositions + sibling->myPositions < 2*logn) //merge !
			{
				// move bits from sibling to x: (more complicated, needs shift!)
				//cout << "merge with left sibling -----------------------------------------" << endl;
				(*x->block)<<=sibling->myPositions;
				x->myPositions+=sibling->myPositions;

				(*x->block)|=(*sibling->block);

				x->myRank+=sibling->myRank;
				updateCountersOnPathToRoot(x);
				rbDelete(sibling, callUpdateCountersOnPathToRoot);
				delete sibling;
			}
			

		
		} // end else if
	} // end else if



}

void BVTree::writeTree(char *writefile){
	std::ofstream write(writefile);
	if (!write)
	{
		cerr << "error: Error writing file: " << writefile<< "." << endl;
		exit(EXIT_FAILURE);
	}
	writeTree(write);
}


size_t* BVTree::getBits(){
	BVNode *n = getRoot();
	size_t blockCounter=0;
	size_t len = getPositions();
	
	int W = sizeof(size_t)*8;
	

	size_t bitsLength = len/W + ((len%W==0)?0:1);

	size_t *bits = new size_t[bitsLength];
	
	size_t i=0; //0 .. bitsLength-1
	int j=0; //0 .. W-1
	
	size_t mask_LeftBitSet = 1;
	
	mask_LeftBitSet <<= W-1;
		
	n=(BVNode*) treeMinimum(root);
	while (n != getNil()) {
		#ifndef NDEBUG
		if (n->block == 0) {
			cerr << "getBits(): block not found !!!" << endl;
			exit(EXIT_FAILURE);
		}
		#endif
		for (blockCounter=0; blockCounter < n->myPositions; blockCounter++) {
			if (j == W) {
				i++;
				j=0;
			}
			bits[i] >>= 1; // set 0

			if ((*n->block)[blockCounter]) bits[i] |= mask_LeftBitSet; // set 1
			j++;	
		}
		n=(BVNode*) treeSuccessor(n);
	}
	
	if (i != bitsLength-1) {
		cout << "last index is wrong" << endl;
		exit(EXIT_FAILURE);
	}
	while (j < W) {
		bits[i] >>= 1;
		j++;
	}
	
	return bits;	
}

void BVTree::writeTree(ostream& stream){
	BVNode *x;
	char c=0;
	size_t cCounter=0;
	size_t blockCounter=0;
	
	x = (BVNode*) treeMinimum(getRoot());

	while (x != getNil()){
		blockCounter=x->myPositions;
		while (blockCounter > 0)
		{
			while (blockCounter > 0 && cCounter < 8 )
			{
				c<<=1;
				if ((*x->block)[x->myPositions - blockCounter]) c++;
				cCounter++;
				blockCounter--;
			}
			if (cCounter == 8) 
			{
				stream << c;
				cCounter = 0;
			}
		}


		x = (BVNode*) treeSuccessor(x);
	}

	if (cCounter != 0) {
		while (cCounter < 8 )
		{
			c<<=1; // fill with zeros.
			cCounter++;
		}
		stream << c;
	}
}


void BVTree::appendBit(bool bit){
	size_t pos = 1;
	if (root != getNil()) pos= getRoot()->subTreePositions + 1;
	
	insertBit(bit, pos);
}

void BVTree::insertBit(bool bit, size_t i){
	if (getRoot() == getNil())
	{
		BVNode *newNode;
		newNode = new BVNode(getNil());
		newNode->color=rbtree::BLACK;
		bitset<2*logn> *newBlock;
		newBlock=new bitset<2*logn>;
		newNode->block = newBlock;
		setRoot(newNode);
		
	}

	BVNode* x = getRoot();
	
	#ifndef NDEBUG
	if (i > x->subTreePositions+1) 
	{
		printNode(x);
		cerr << "error: insertBit: position " << i <<" in block not available, only " << x->myPositions <<" positions.\n"; //shouldn't happen
		exit(EXIT_FAILURE);
	}
	if (i < 1) 
	{
		cerr << "error: insertBit: position " << i <<" is not valid.\n"; //shouldn't happen
		exit(EXIT_FAILURE);
	}
	#endif

	size_t lsP=0;
	
	
	while (true) 
	{
		// update subTree-counters
		(x->subTreePositions)++;
		if (bit) (x->subTreeRank)++;
		
		lsP = x->getLeft()->subTreePositions;
		
		if (lsP >= i) 
		{//	cout << "A" << endl;
			x=x->getLeft();
		}
		else if (lsP+x->myPositions >= i-1){ //-1 to append to last position in current node
			i-=lsP;
			break;
		}
		else{
			i-=(lsP+x->myPositions);
			x=x->getRight();
		}
	}

	// now put the bit into the block x:
	bitset<2*logn> mask;
	
	if (x->myPositions > 0) 
	{
		if (i != 1) 
		{
			mask.set();
			mask>>=2*logn - i +1; // make 1's at all positions < i
			mask &= *(x->block);      // store bits in front of new bit	
			(*(x->block)>>=i-1)<<=i;  // shift bits by one
			(*x->block)|=mask;        // restore bits in front of new bit
		}
		else 
		  (*x->block)<<=1;
	}
	
	(*x->block)[i-1]=bit;     // set new bit

	// update local pointers
	(x->myPositions)++;      //update p (leaf)
	if (bit) (x->myRank)++;  //update r (leaf)

	#ifndef NDEBUG	
	if ((int)x->myPositions > 2*logn) 
	{
		cerr << "error: positions in block already too many.\n"; //shouldn't happen
		exit(EXIT_FAILURE);
	}
	#endif

	// split and rotate, if necessary
	if ((int)x->myPositions == 2*logn)
	{
		//cout << "split !-------------" << endl;

		// new node:		
		BVNode *newNode;
		newNode = new BVNode(getNil());

		//find place for new node:
		BVNode *y;// some parent of new node
		if (x->getRight() != getNil()) {
			y = (BVNode*) this->treeMinimum(x->getRight());
			y->setLeft(newNode);
		} else {
			y=x;
			y->setRight(newNode);
		}
		newNode->setParent(y);
		
		//new block:
		bitset<2*logn> *newBlock;
		newBlock=new bitset<2*logn>;
		*newBlock=*x->block>>logn; //copy of bits into the new block
		
		mask.set();
		mask>>=logn;
		*x->block &= mask; //delete bits that already have been copied to the new block		
		
		newNode->block=newBlock;
		newNode->myRank=(*newNode->block).count();
		newNode->myPositions=logn;
		newNode->color=rbtree::RED;

		//update old node x:
		x->myRank-=newNode->myRank;     //=(*x->block).count();
		x->myPositions=logn;

		updateCountersOnPathToRoot(newNode);  // meets x
		rbInsertFixup(newNode, callUpdateCounters);
		
	} // end if
	

}

void BVTree::checkSubTree(BVNode *n){
	size_t lP = 0;
	size_t lR = 0;
	size_t rP = 0;
	size_t rR = 0;

	if (n->getLeft()!=getNil()) {
		lP = n->getLeft()->subTreePositions;
		lR = n->getLeft()->subTreeRank;
		if (n->getLeft()->getParent() != n) {
			cout << "au"<< endl;
			exit(1);
		}
		
	}

	if (n->getRight()!=getNil()) {
		rP = n->getRight()->subTreePositions;
		rR = n->getRight()->subTreeRank;
		if (n->getRight()->getParent() != n) {
			cout << "au"<< endl;
			exit(1);
		}		
	}

	if  ( (n->subTreePositions != (n->myPositions + lP + rP)) ||
	      (n->subTreeRank != (n->myRank + lR + rR)) ){
		cout << "checkSubTree: error" <<endl;
		cout << "lP: " << lP << endl;
		cout << "lR: " << lR << endl;
		cout << "rP: " << rP << endl;
		cout << "rR: " << rR << endl;
		cout << "n->myPositions + lP + rP: " << n->myPositions + lP + rP << endl;
		cout << "n->myRank + lR + rR: " << n->myRank + lR + rR << endl;
		printNode(n);
		printNode(n->getLeft());
		printNode(n->getRight());
		exit(1);
	} 	

	if (n->getLeft()!=getNil()) checkSubTree(n->getLeft());
	if (n->getRight()!=getNil()) checkSubTree(n->getRight());
}

void callUpdateCounters(RBNode *n, RBTree *T){
	((BVTree*)T)->updateCounters((BVNode*) n);
}

void callUpdateCountersOnPathToRoot(RBNode *n, RBTree *T){
	((BVTree*)T)->updateCountersOnPathToRoot((BVNode*) n);
}

} // namespace

