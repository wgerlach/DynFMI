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


 //TODO check star color of the nodes pos and handle!
 
 using namespace std;

#ifndef GUARD_PosTree
#define GUARD_PosTree

#ifndef uchar
#define uchar unsigned char
#endif
#ifndef size_t
#define size_t unsigned long
#endif

#include <stdio.h>
#include <stdlib.h>

  
#include "rbtree.h"
#include "handle.h"


namespace handletree{
class HandleTree;
class HandleNode;
}


namespace postree {

//typedef rbtree::RBNode RBNode;
//typedef rbtree::RBTree RBTree;
//typedef rbtree::BLACK BLACK;
//typedef rbtree::RED RED;
using namespace rbtree;



typedef handletree::HandleTree HandleTree;
typedef handletree::HandleNode HandleNode;


//class Handle;
//class HandleNode;

void callPosUpdateCounters(RBNode *n, RBTree *T);
void callPosUpdateCountersOnPathToRoot(RBNode *n, RBTree *T);


 
class PosNode : public RBNode {
	public:
	
	HandleNode *handleNode;
	size_t subtreeSize;
	size_t textSize; // size including endmarker!
	//size_t sumTextSize; // sum of textlength of all texts located in this subtree;
	
	PosNode() // NIL
		: RBNode(this),subtreeSize(0),textSize(0) { 
	}

	PosNode(PosNode *n, size_t textSize)
		: RBNode(n),subtreeSize(1),textSize(textSize) { 
	}
	
		
	PosNode* getParent(){
		return ((PosNode*) ((RBNode*) this)->parent);
	}

	PosNode* getLeft(){
		return ((PosNode*) ((RBNode*) this)->left);
	}

	PosNode* getRight(){
		return ((PosNode*) ((RBNode*) this)->right);
	}

	void setParent(PosNode* n){
		((RBNode*) this)->parent=(RBNode*)n;
	}

	void setLeft(PosNode* n){
		((RBNode*) this)->left=(RBNode*)n;
	}

	void setRight(PosNode* n){
		((RBNode*) this)->right=(RBNode*)n;
	}	
};

class PosTree: public RBTree
{
	public:
	
	HandleTree* handle;
	
	size_t textNumber;
	size_t matchPosition;

	size_t sampleInterval;

	PosTree(size_t sampleInterval){
		setNil(new PosNode());
		setRoot(getNil());
		this->sampleInterval=sampleInterval;
	}

	void setRoot(PosNode* n){
		((RBTree*) this)->root=(RBNode*)n;
	}
	
	PosNode* getRoot(){
		return ((PosNode*) ((RBTree*) this)->root);
	}

	void setNil(PosNode* n){
		((RBTree*) this)->nil=(RBNode*)n;
	}

	PosNode* getNil(){
		return ((PosNode*) ((RBTree*) this)->nil);
	}

	size_t getSize(){
		return (getRoot()!=getNil())?getRoot()->subtreeSize:0;
	}


	
	size_t getPos(PosNode *n);


	PosNode* getPosNode(size_t text);
	size_t getTextSize(size_t pos);
	void deleteText(size_t pos);
	void deletePosNode(PosNode *n);
	size_t appendText(size_t textSize);

	void updateCountersOnPathToRoot(PosNode *n);
	void updateCounters(PosNode *n);
	
};

} // namespace


#endif
