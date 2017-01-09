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


#include "handle.h"


namespace handletree {

size_t HandleTree::getPos(size_t key){
	HandleNode *n = getKey(key);
	if (n==0) return 0;

	PosNode *p = n->posNode;

	#ifndef NDEBUG
	if (p==pos->getNil()) {
		cerr << "error: HandleTree::getPos: no PosNode found." << endl;
		exit(EXIT_FAILURE);
	}
	#endif
	return (this->pos->getPos(p));
}

size_t* HandleTree::getKeys(){
	size_t num = getSize();
	size_t i = 0;
	
	if (this->root == getNil()) return 0;
	
	size_t *keys= new size_t[num];
	
	
	HandleNode *n = (HandleNode*) treeMinimum(root);
	#ifndef NDEBUG
	HandleNode *oldn;
	#endif
	while (n != getNil()) {
		
		//printf("node: %d", n);
		keys[i] = n->key;
		i++;
		#ifndef NDEBUG
		oldn = n;
		#endif
		n = (HandleNode*) treeSuccessor((RBNode*) n);
		#ifndef NDEBUG
		if (n == oldn) {
				cerr << "old!" << endl;
				exit(EXIT_FAILURE);
				}
		#endif
	}
	#ifndef NDEBUG
	if (num != i) {
		cerr << "num: " << num << endl;
		cerr << "i: " << i << endl;
		cerr << "error: HandleTree::getKeys: key number was not correct." << endl;
		exit(EXIT_FAILURE);
	}
	#endif
	return keys;
}

void HandleTree::deleteKey(size_t key){

	HandleNode *n = getKey(key);
	this->pos->deletePosNode(n->posNode);
	
	rbDelete( n, callHandleUpdateCountersOnPathToRoot);
	delete n;
}




HandleNode* HandleTree::getKey(size_t key){
	HandleNode *n= getRoot();
	while (n != getNil()) {
		if (n->key == key) return n;
		if (n->getLeft() != getNil()) {
			if (key <= n->getLeft()->maxkey) n=n->getLeft();
			else n=n->getRight();
		}
		else n=n->getRight();
	}
	
			
	return 0;
}

void HandleTree::updateCountersOnPathToRoot(HandleNode *n) {
	while (n != getNil()) {
		updateCounters(n);
		
		n = n->getParent();
	}
}


void HandleTree::updateCounters(HandleNode *n) {
	#ifndef NDEBUG
	if (n == getNil()) {
		cerr << "error: HandleTree::updateCounters" << endl;
		exit(EXIT_FAILURE);
	}
	#endif
	n->subtreeSize=1;
	
	
	n->subtreeSize += n->getLeft()->subtreeSize;
	
	
	if ( n->getRight() != getNil()) {
		n->subtreeSize += n->getRight()->subtreeSize;
		n->maxkey=n->getRight()->maxkey;
	} else n->maxkey=n->key;
	
}

		
HandleNode* HandleTree::getNewKey(){
	HandleNode *n= getRoot();
	
	if (n == getNil()) {
		//tree is empty
		HandleNode *newLeaf= new HandleNode(getNil(),1); // 1=smallest key, 0 = error
		setRoot(newLeaf);
		((RBNode*)newLeaf)->color=BLACK;
		
		return newLeaf;
	}
	
	if (n->maxkey == n->subtreeSize) {
		// tree is full
		HandleNode *last = (HandleNode*) treeMaximum(n);
		
		
		HandleNode *newLeaf= new HandleNode(getNil(), n->maxkey+1);
		
		last->setRight(newLeaf);
		newLeaf->setParent(last);
		
		if (newLeaf->getParent() != getNil()) updateCountersOnPathToRoot(newLeaf->getParent());
		
		rbInsertFixup(newLeaf, callHandleUpdateCounters);
		
		return newLeaf;
	}
	
	HandleNode *newNode;
	size_t smallerkeys = 0;
	size_t lmax; //getLeft()->maxkey
	size_t lsub; //n->getLeft()->subtreeSize
	while (true) {
		#ifndef NDEBUG
		cout << "search first free key" << endl;
		#endif
		// search first free key
		
		lmax = n->getLeft()->maxkey;
		lsub = n->getLeft()->subtreeSize;
		

		if (lmax == 0) { // no left child
			if (smallerkeys+1 < n->key) { // check if it is free
				newNode= new HandleNode(getNil(), smallerkeys+1);
				newNode->setParent(n);
				n->setLeft(newNode);
				updateCountersOnPathToRoot(n);
				rbInsertFixup(newNode, callHandleUpdateCounters);
				return newNode;
			}
		} else { //left child exists
			if ( lmax > (lsub + smallerkeys) ) { 
				// free key at left subtree
				n=n->getLeft();
				continue;
			} else if (lmax + 1 < n->key) { // full left subtree, check if it is free inbetween
				// insert to predecessor
				n=(HandleNode*)treePredecessor(n);
				//found place ->  :predeccessor new key = lmax + 1
				newNode= new HandleNode(getNil(), lmax+1);
				newNode->setParent(n);
				n->setRight(newNode);
				updateCountersOnPathToRoot(n);
				rbInsertFixup(newNode, callHandleUpdateCounters);
				return newNode;
	
			}
		}
		
		smallerkeys += 1+lsub;

		if (n->getRight() == getNil()) { // no right child, must be free
			newNode= new HandleNode(getNil(), smallerkeys+1);
			newNode->setParent(n);
			n->setRight(newNode);
			updateCountersOnPathToRoot(n);
			rbInsertFixup(newNode, callHandleUpdateCounters);
			return newNode;
		} else { //right child exists
			
			//n=n->getRight();
			if (n->getRight()->maxkey-smallerkeys == n->getRight()->subtreeSize) { // right subTree is full, insert in front
				size_t leftMinKey = n->getRight()->maxkey - n->getRight()->subtreeSize;
				if (leftMinKey -1 >  n->key) { //in front there is space
					//insert new n-key + 1
					newNode= new HandleNode(getNil(), n->key+1);
					n=(HandleNode*)treeSuccessor(n);
					newNode->setParent(n);
					n->setLeft(newNode);
					updateCountersOnPathToRoot(n);
					rbInsertFixup(newNode, callHandleUpdateCounters);
					return newNode;
				} else {
					cerr << "error: HandleTree::getNewKey: no space ?!? " << endl;
					exit(EXIT_FAILURE);
				}
			} else {
				n=n->getRight();
			}


		}
	
		#ifndef NDEBUG	
		if (n==getNil()) {
			cerr << "error: HandleTree::getNewKey: (A) something wrong ! " << endl;
			exit(EXIT_FAILURE);
		}
		#endif
	}
	#ifndef NDEBUG
	cerr << "error: HandleTree::getNewKey: (B) something wrong ! " << endl;
	exit(EXIT_FAILURE);
	#endif
	return 0; // error
}

void callHandleUpdateCounters(RBNode *n, RBTree *T){
	((HandleTree *)T)->updateCounters((HandleNode*) n);
}

void callHandleUpdateCountersOnPathToRoot(RBNode *n, RBTree *T){
	((HandleTree *)T)->updateCountersOnPathToRoot((HandleNode*) n);
}

} // namespace

