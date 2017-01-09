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

// Implementation of the Dynamic Bit Vector with Indels problem
// space: O(n) time: O(log(n))
// papers: V. Maekinen, G. Navarro. Dynamic Entropy-Compressed Sequences and Full-Text
//           Indexes. CPM 2006, Chapter 3.4 Dynamic Structures for Bit Vectors
//   also: W.-L. Chan, W.-K. Hon, and T.-W. Lam. Compressed index for a dynamic collection
//           of texts. In Proc. CPM04, LNCS 3109, pages 445-456, 2004 

#ifndef GUARD_BVTree
#define GUARD_BVTree

#include <iostream>
#include <bitset>
#include <cstdlib>
//#include <map>
//#include <stack>
#include <cmath>
#include <fstream>
//#include <vector>
#include <cstdio>

#include "rbtree.h"


#ifndef uchar
#define uchar unsigned char
#endif
#ifndef size_t
#define size_t unsigned long
#endif


#ifndef LOGN
#define LOGN 64
#endif

const int logn = LOGN;

//upperBound = 2 * logn;
//lowerBound = logn / 2;


class BVNode;
class BVTree;


namespace bvtree {

typedef rbtree::RBNode RBNode;
typedef rbtree::RBTree RBTree;


void callUpdateCounters(RBNode *n, RBTree *T);
void callUpdateCountersOnPathToRoot(RBNode *n, RBTree *T);

class BVNode : public RBNode
{
	public:
	size_t myPositions; // 4*4 bytes = 16 bytes
	size_t myRank;
	size_t subTreePositions; //number of positions stored in the subtree rooted at this node
	size_t subTreeRank;      //number of bits set in the subtree rooted at this node

	bitset<2*logn> *block; // 4 bytes
	

	BVNode()
		: RBNode(this), myPositions(0), myRank(0), subTreePositions(0), subTreeRank(0), block(0) {
	}


	BVNode(BVNode* n)
		: RBNode(n), myPositions(0), myRank(0), subTreePositions(0), subTreeRank(0), block(0) {
	}

	~BVNode(){
		delete block;
	}

		
	BVNode* getParent(){
		return ((BVNode*) ((RBNode*) this)->parent);
	}

	BVNode* getLeft(){
		return ((BVNode*) ((RBNode*) this)->left);
	}

	BVNode* getRight(){
		return ((BVNode*) ((RBNode*) this)->right);
	}

	void setParent(BVNode* n){
		((RBNode*) this)->parent=(RBNode*)n;
	}

	void setLeft(BVNode* n){
		((RBNode*) this)->left=(RBNode*)n;
	}

	void setRight(BVNode* n){
		((RBNode*) this)->right=(RBNode*)n;
	}
		
};


class BVTree : public RBTree{
public:

  //Constructors
  BVTree(){
	
	setNil(new BVNode());
	setRoot(getNil());

  	tempnil = getNil();

  	tempbit = new bitset<2*logn>;
  }

  //Destructor:
  ~BVTree();

  bool operator[](size_t);


  // inserts bit at position i, countings begins with 1:
  void appendBit(bool bit);
  void insertBit(bool bit, size_t i);
  void deleteBit(size_t i);

  size_t rank0(size_t i);
  size_t rank1(size_t i);
  size_t rank(bool b, size_t i){return b?rank1(i):rank0(i);}
  
  size_t select0(size_t i);
  size_t select1(size_t i);
  size_t select(bool b, size_t i){return b?select1(i):select0(i);}

  void setRoot(BVNode* n){
	((RBTree*) this)->root=(RBNode*)n;
  }
	
  BVNode* getRoot(){
	return ((BVNode*) ((RBTree*) this)->root);
  }

  void setNil(BVNode* n){
	tempnil = n;
	((RBTree*) this)->nil=(RBNode*)n;
  }

  BVNode* getNil(){
  	return tempnil;
	
  }

  // write bits back into a stream:  
  size_t* getBits();
  void writeTree(char *writefile); 
  void writeTree(std::ostream& stream); //e.g. stream = cout

  int getTreeMaxDepth();
  int getTreeMinDepth();
  size_t getPositions();
  size_t getRank();
  
  void iterateReset();
  bool iterateGetBit();
  bool iterateNext();
  size_t iterateGetRank(bool bit);
  
  bool getLastBitDeleted(){return lastBitDeleted;}
  size_t getLastRank(){return lastRank;}
  
  void checkSubTree(BVNode *n);

  void updateCounters(BVNode *n);
  void updateCountersOnPathToRoot(BVNode *walk);
  
  //debug:
  void printNode(size_t i);

protected:

  size_t iterate;
  size_t iterateLocal;
  size_t iterateRank;
  BVNode *iterateNode;

  BVNode *tempnil;

  bool lastBitDeleted;
  size_t lastRank;
  

  bitset<2*logn> *tempbit;

  // content of BVNode, for debugging:
  void printNode(BVNode *n);
  
  // other operations:
  size_t getLocalRank(BVNode* n, size_t position);
  size_t getLocalSelect0(BVNode* n, size_t query);
  size_t getLocalSelect1(BVNode* n, size_t query);
  
  void deleteNode(BVNode *n);
  void deleteLeaf(BVNode *leaf);
};

} // namespace

#endif

