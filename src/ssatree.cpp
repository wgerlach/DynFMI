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

#include "ssatree.h"

namespace sampledsatree {

void SampledSATree::getSample(size_t i){
	SampledSANode *n= getRoot();
	size_t rest=0;
	size_t left_subtreeSize;


	while (true) {
		#ifndef NDEBUG
		if (n == getNil()) {
				cout << "getSample: nil!" << endl;
				exit(EXIT_FAILURE);
			}
		#endif

		left_subtreeSize = n->getLeft()->subtreeSize;
	
		if (i <= left_subtreeSize + rest) {
			n=n->getLeft();
		} else if (i > left_subtreeSize + rest + 1) {
			rest += left_subtreeSize + 1;
			n= n->getRight();
		} else {
			// match
			this->handle=n->handle;
			this->sampledSAValue = n->sampledSAValue;
			return;
		}
	}

}

void SampledSATree::insertSample(size_t sample, size_t handle, size_t i){
	SampledSANode *n= getRoot();
	
	
	if (n == getNil()) {
		SampledSANode *newLeaf= new SampledSANode(getNil(), handle, sample); // 1=smallest key, 0 = error
		setRoot(newLeaf);
		((RBNode*)newLeaf)->color=BLACK;
		return;
	}
	

	if (i > n->subtreeSize ) {
		//append
		n= (SampledSANode*) treeMaximum(n);
		SampledSANode *newLeaf= new SampledSANode(getNil(), handle, sample);
		n->setRight(newLeaf);
		newLeaf->setParent(n);
		updateCountersOnPathToRoot(newLeaf);
		#ifndef NDEBUG
		if (newLeaf == getNil() || newLeaf == 0) {
			cout << "insertSample, newLeaf: nil!" << endl;
			exit(EXIT_FAILURE);
			}
		#endif
		rbInsertFixup(newLeaf, callSampledSATreeUpdateCounters);
		return;
	}

	size_t rest=0;
	size_t left;


	while (true) {
		#ifndef NDEBUG
		if (n == getNil()) {
				cout << "insertSample: nil!" << endl;
				exit(EXIT_FAILURE);
			}
		#endif

		left = n->getLeft()->subtreeSize;
		
	
		if (i <= left + rest) {
			n=n->getLeft();
		} else if (i == left + rest + 1) {
			
			// match
			SampledSANode *newLeaf= new SampledSANode(getNil(), handle, sample);
			if (n->getLeft() == getNil()) {
				n->setLeft(newLeaf);
				
			} else {
				n=(SampledSANode*) treeMaximum(n->getLeft());
				n->setRight(newLeaf);

			}
			newLeaf->setParent(n);
			#ifndef NDEBUG
			if (newLeaf == getNil() || newLeaf == 0) {
				cout << "insertSample, newLeaf_2: nil!" << endl;
				exit(EXIT_FAILURE);
			}
			#endif
			
			updateCountersOnPathToRoot(newLeaf);
			rbInsertFixup(newLeaf, callSampledSATreeUpdateCounters);
			return;
		} else {
			rest+=left+1;
			n=n->getRight();
		}
	}
	
	
}

void SampledSATree::deleteSample(size_t i){
	SampledSANode *n= getRoot();
	
	

	size_t rest=0;
	size_t left_subtreeSize;


	while (true) {
		left_subtreeSize = n->getLeft()->subtreeSize;
		
	
		if (i < left_subtreeSize + rest) {
			n=n->getLeft();
		} else if (i > left_subtreeSize + rest) {
			rest += left_subtreeSize + 1;
			n= n->getRight();
		} else {
			// match
			rbDelete(n, callSampledSATreeUpdateCountersOnPathToRoot);
			delete n;
		}
	}
	
	
}


void SampledSATree::printSubTree(SampledSANode *n){
	cout << n->handle << "/" << n->sampledSAValue << " [";
	if (n->getLeft()==getNil()) cout << "N";
		else printSubTree(n->getLeft());
	cout << ",";
	if (n->getRight()==getNil()) cout << "N";
		else printSubTree(n->getRight());
	cout << "]";
}


void SampledSATree::updateCountersOnPathToRoot(SampledSANode *n) {
	while (n != getNil()) {
		updateCounters(n);
		
		n = n->getParent();
	}
}


void SampledSATree::updateCounters(SampledSANode *n) {
	
	if (n == getNil()) {
		#ifndef NDEBUG
		cerr << "error: SampledSATree::updateCounters" << endl;
		exit(EXIT_FAILURE);
		#endif	
		return;
		
	}
	

	n->subtreeSize=1;
	
	n->subtreeSize += n->getLeft()->subtreeSize; 	
	
	n->subtreeSize += n->getRight()->subtreeSize;
}

void callSampledSATreeUpdateCounters(RBNode *n, RBTree *T){
	#ifndef NDEBUG
	if (n == T->nil) {
		cerr << "error: callSampledSATreeUpdateCounters n=NIL" << endl;
		exit(EXIT_FAILURE);	
	}
	#endif	
	((SampledSATree *)T)->updateCounters((SampledSANode*) n);
}

void callSampledSATreeUpdateCountersOnPathToRoot(RBNode *n, RBTree *T){
	#ifndef NDEBUG
	if (n == T->nil) {
		cerr << "error: callSampledSATreeUpdateCountersOnPathToRoot n=NIL" << endl;
	exit(EXIT_FAILURE);	
	}
	#endif		
	((SampledSATree *)T)->updateCountersOnPathToRoot((SampledSANode*) n);
}


} //namespace




