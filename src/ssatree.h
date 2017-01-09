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


 using namespace std;

#ifndef GUARD_SampledSATree
#define GUARD_SampledSATree
 
#ifndef uchar
#define uchar unsigned char
#endif
#ifndef size_t
#define size_t unsigned long
#endif

#include <stdio.h>
#include <stdlib.h>


  
#include "rbtree.h"


namespace sampledsatree {


//typedef rbtree::RBNode RBNode;
//typedef rbtree::RBTree RBTree;
//typedef rbtree::RBNodecolor RBNodecolor;

using namespace rbtree;

class SampledSATree;
class SampledSANode;

void callSampledSATreeUpdateCounters(RBNode *n, RBTree *T);
void callSampledSATreeUpdateCountersOnPathToRoot(RBNode *n, RBTree *T);

 
class SampledSANode : public RBNode {
	public:
	size_t sampledSAValue; //maximum for internal nodes
	size_t subtreeSize;
	size_t handle;

	SampledSANode() 
		: RBNode(this), sampledSAValue(0), subtreeSize(0), handle(0){
	}	


	SampledSANode(SampledSANode *n, size_t handle, size_t sampledSAValue) 
		: RBNode(n), sampledSAValue(sampledSAValue), subtreeSize(1), handle(handle){
	}	
		
	SampledSANode* getParent(){
		return ((SampledSANode*) ((RBNode*) this)->parent);
	}

	SampledSANode* getLeft(){
		return ((SampledSANode*) ((RBNode*) this)->left);
	}

	SampledSANode* getRight(){
		return ((SampledSANode*) ((RBNode*) this)->right);
	}

	void setParent(SampledSANode* n){
		((RBNode*) this)->parent=(RBNode*)n;
	}

	void setLeft(SampledSANode* n){
		((RBNode*) this)->left=(RBNode*)n;
	}

	void setRight(SampledSANode* n){
		((RBNode*) this)->right=(RBNode*)n;
	}
	
};

class SampledSATree : public RBTree{
	public:

	//results:
	size_t sampledSAValue;
	size_t handle;

	SampledSATree(){
		setNil(new SampledSANode());
		setRoot(getNil());
	}

	void setRoot(SampledSANode* n){
		((RBTree*) this)->root=(RBNode*)n;
	}
	
	SampledSANode* getRoot(){
		return ((SampledSANode*) ((RBTree*) this)->root);
	}

	void setNil(SampledSANode* n){
		((RBTree*) this)->nil=(RBNode*)n;
	}

	SampledSANode* getNil(){
		return ((SampledSANode*) ((RBTree*) this)->nil);
	}


	size_t getSize(){
		return (getRoot()!=getNil())?getRoot()->subtreeSize:0;
	}
	

	void getSample(size_t i); // sets sampledSAValue and handle
	void insertSample(size_t sample, size_t handle, size_t i);
	void deleteSample(size_t i);


	void updateCountersOnPathToRoot(SampledSANode *n);
	void updateCounters(SampledSANode *n);
	void printSubTree(SampledSANode *n);
	
};

} // namespace

#endif
