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

#include <iostream>

#include "pos.h"

namespace postree {

size_t PosTree::getPos(PosNode *n){
	size_t position=0;
	
	if (n->getLeft() != getNil()) position += n->getLeft()->subtreeSize;
	
	while (n->getParent() != getNil()) {
		if (n == n->getParent()->getRight()) {
			// rightchild
			if (n->getParent()->getLeft() != getNil()) position += n->getParent()->getLeft()->subtreeSize;
			position++;
		} 
		
		n=n->getParent();
	}
	
	
	return position+1;
}

size_t PosTree::getTextSize(size_t pos){
	PosNode *n= getPosNode(pos);
	
	return n->textSize;
}

void PosTree::deleteText(size_t pos){
	PosNode *n= getPosNode(pos);
	
	// current node matches.
	rbDelete(n, callPosUpdateCountersOnPathToRoot);
	delete n;
}

void PosTree::deletePosNode(PosNode *n){
	rbDelete(n, callPosUpdateCountersOnPathToRoot);
	delete n;
}

PosNode* PosTree::getPosNode(size_t pos){
	// smallest pos is 1 !
	pos--;
	PosNode *n= getRoot();
	size_t leftTree=0;
	size_t leftSubTreeSize;
	while (true) {
		leftSubTreeSize = n->getLeft()->subtreeSize;
		
		if (pos == leftTree + leftSubTreeSize ) {
			// current node matches.
			return n;
		} else if (pos < leftTree + leftSubTreeSize){
			n=n->getLeft();
			}
		else {
			leftTree += leftSubTreeSize +1;
			n=n->getRight();
			}
	}
	
	
	cerr << "error: PosTree::getPosNode: text POS " << pos << " not found!" << endl;
	exit(EXIT_FAILURE);
			
	return  getNil();
}


size_t PosTree::appendText(size_t textSize){
	PosNode *n= getRoot();
	
	if (n == getNil()) {
		PosNode *newLeaf= new PosNode(getNil(), textSize);

		setRoot(newLeaf);
		((RBNode*)newLeaf)->color = BLACK;
		
		HandleNode* hn = handle->getNewKey();
	
		// connect handleNode and posNode:
		hn->posNode = newLeaf;
		newLeaf->handleNode = hn;
	
		return hn->key;
	}
	
	while (n->getRight() != getNil()) {
		n=n->getRight();
	}

	// new right child !
	PosNode *newLeaf= new PosNode(getNil(), textSize);


	n->setRight(newLeaf);
	newLeaf->setParent(n);
	
	updateCountersOnPathToRoot(newLeaf);

	rbInsertFixup(newLeaf, callPosUpdateCounters);

	HandleNode* hn = handle->getNewKey();
	
	// connect handleNode and posNode:
	hn->posNode = newLeaf;
	newLeaf->handleNode = hn;
	
	return hn->key;
}




void PosTree::updateCountersOnPathToRoot(PosNode *n) {
	while (n != getNil()) {
		updateCounters(n);
		
		n = n->getParent();
	}
}

void PosTree::updateCounters(PosNode *n) {
	n->subtreeSize=1;
	//n->sumTextSize = n->textSize;
	
	n->subtreeSize += n->getLeft()->subtreeSize;
	//n->sumTextSize += n->getLeft()->sumTextSize;
	
	n->subtreeSize += n->getRight()->subtreeSize;
	//n->sumTextSize += n->getRight()->sumTextSize;
}



void callPosUpdateCounters(RBNode *n, RBTree *T){
	((PosTree *)T)->updateCounters((PosNode*) n);
}

void callPosUpdateCountersOnPathToRoot(RBNode *n, RBTree *T){
	((PosTree *)T)->updateCountersOnPathToRoot((PosNode*) n);
}

} // namespace
