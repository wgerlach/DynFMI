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

#ifndef GUARD_HandleTree
#define GUARD_HandleTree
 
#ifndef uchar
#define uchar unsigned char
#endif
#ifndef size_t
#define size_t unsigned long
#endif
  
#include "rbtree.h"
#include "pos.h"


namespace postree{
class PosTree;
class PosNode;
}

namespace handletree {

//typedef rbtree::RBNode RBNode;
//typedef rbtree::RBTree RBTree;
//typedef rbtree::RBNodecolor RBNodecolor;
//typedef rbtree::BLACK BLACK;
//typedef rbtree::RED RED;

using namespace rbtree;

typedef postree::PosTree PosTree;
typedef postree::PosNode PosNode;

void callHandleUpdateCounters(RBNode *n, RBTree *T);
void callHandleUpdateCountersOnPathToRoot(RBNode *n, RBTree *T);

 
class HandleNode : public RBNode {
	public:
	size_t key; //maximum for internal nodes
	size_t subtreeSize;
	size_t maxkey;
	PosNode *posNode;

	HandleNode() 
		: RBNode(this), key(0), subtreeSize(0), maxkey(0){
	}


	HandleNode(HandleNode *n, size_t key) 
		: RBNode(n), key(key), subtreeSize(1), maxkey(key){
	}
	
	
	HandleNode* getParent(){
		return ((HandleNode*) ((RBNode*) this)->parent);
	}

	HandleNode* getLeft(){
		return ((HandleNode*) ((RBNode*) this)->left);
	}

	HandleNode* getRight(){
		return ((HandleNode*) ((RBNode*) this)->right);
	}

	void setParent(HandleNode* n){
		((RBNode*) this)->parent=(RBNode*)n;
	}

	void setLeft(HandleNode* n){
		((RBNode*) this)->left=(RBNode*)n;
	}

	void setRight(HandleNode* n){
		((RBNode*) this)->right=(RBNode*)n;
	}
	
};

class HandleTree : public RBTree{
	public:
	PosTree *pos;

	HandleTree(){
		setNil(new HandleNode());
		setRoot(getNil());
	}
	

	void setRoot(HandleNode* n){
		((RBTree*) this)->root=(RBNode*)n;
	}
	
	HandleNode* getRoot(){
		return ((HandleNode*) ((RBTree*) this)->root);
	}

	void setNil(HandleNode* n){
		((RBTree*) this)->nil=(RBNode*)n;
	}

	HandleNode* getNil(){
		return ((HandleNode*) ((RBTree*) this)->nil);
	}


	size_t getSize(){
		return (getRoot()!=getNil())?getRoot()->subtreeSize:0;
	}
	
	size_t* getKeys();
	
	size_t getPos(size_t key);
	
	HandleNode* getKey(size_t key);
	void deleteKey(size_t key);
	void updateCountersOnPathToRoot(HandleNode *n);
	void updateCounters(HandleNode *n);
		
	HandleNode* getNewKey();
	
	
	
};

} // namespace

#endif
