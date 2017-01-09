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
#include <cstdlib>
#include "rbtree.h"


using namespace std;

namespace rbtree {

void RBTree::checkTree(){
	cout << "start" << endl;
	this->root->countBlack(0);
	cout << "stop" << endl;
}

void RBTree::leftRotate(RBNode* x, void (*updateNode)(RBNode* n, RBTree *T)){
	RBNode* y = x->right;
	x->right=y->left;
	if (y->left != nil) y->left->parent=x;
	y->parent=x->parent;
	if (x->parent==this->nil)
		root=y;
	else {
		if (x == x->parent->left)
			x->parent->left=y;
		else
		x->parent->right=y;
		}
	y->left=x;
	x->parent=y;

	// update counters of x and y !
	if (updateNode != 0) {
		updateNode(x, this);
		updateNode(y, this);
		}
}

void RBTree::rightRotate(RBNode* x, void (*updateNode)(RBNode* n, RBTree *T)){
	RBNode* y = x->left;
	x->left=y->right;
	if (y->right != nil) y->right->parent=x;
	y->parent=x->parent;


	if (x->parent==nil)
		root=y;
	else {
		if (x == x->parent->right) {
			x->parent->right=y;
		} else {
			x->parent->left=y;
		} 
	}

	y->right=x;
	x->parent=y;
	if (updateNode != 0) {
		updateNode(x, this);
		updateNode(y, this);
		}
}

void RBTree::rbInsertFixup(RBNode* z, void (*updateNode)(RBNode* n, RBTree *T)){
	RBNode *y;
	
	if (z->parent==nil) return;
	while (z->parent->color == RED) {
		if (z->parent == z->parent->parent->left)
		{
			y = z->parent->parent->right;
			
			if (y->color==RED)
			{
				z->parent->color=BLACK;        // Case 1a
				y->color=BLACK;                // Case 1a
				z->parent->parent->color=RED;  // Case 1a
				z=z->parent->parent;           // Case 1a
			}
			else
			{
				if (z==z->parent->right)
				{
					z=z->parent;      // Case 2a
					leftRotate(z, updateNode);  // Case 2a
				}
				z->parent->color=BLACK;             // Case 3a
				z->parent->parent->color=RED;       // Case 3a
				rightRotate(z->parent->parent, updateNode);  // Case 3a
			}
		}
		else
		{

			y = z->parent->parent->left;

			if (y->color==RED)
			{
				z->parent->color=BLACK;        // Case 1b
				y->color=BLACK;                // Case 1b
				z->parent->parent->color=RED;  // Case 1b
				z=z->parent->parent;           // Case 1b
			}
			else
			{
				if (z==z->parent->left)
				{
					z=z->parent;      // Case 2b
					rightRotate(z, updateNode);  // Case 2b
				}
				z->parent->color=BLACK;             // Case 3b
				z->parent->parent->color=RED;       // Case 3b
				leftRotate(z->parent->parent, updateNode);  // Case 3b
			}
		}
		if (z->parent==nil) return;
	} //end of while

	root->color=BLACK;
}


void RBTree::rbDeleteFixup(RBNode *x, void (*updateNode)(RBNode* n, RBTree *T)){
	RBNode *w;

	while ((x != root) && (x->color == BLACK))
	{
		if (x == x->parent->left)
		{
			w=x->parent->right;
			
			if (w->color == RED)
			{
				w->color=BLACK;           // Case 1a
				x->parent->color=RED;     // Case 1a
				leftRotate(x->parent, updateNode);  // Case 1a
				w=x->parent->right;       // Case 1a
			}
			if ((w->left->color == BLACK) && (w->right->color == BLACK))
			{
				w->color=RED;  // Case 2a
				x=x->parent;   // Case 2a
			}
			else
			{
				if (w->right->color == BLACK)
				{
					w->left->color=BLACK;  // Case 3a
					w->color=RED;          // Case 3a
					rightRotate(w, updateNode);      // Case 3a
					w=x->parent->right;    // Case 3a
				}
				w->color=x->parent->color;  // Case 4a
				x->parent->color=BLACK;     // Case 4a
				w->right->color=BLACK;      // Case 4a
				leftRotate(x->parent, updateNode);    // Case 4a
				x=root;                  // Case 4a
			}
		}
		else
		{
			w=x->parent->left;

			if (w->color == RED)
			{
				w->color=BLACK;           // Case 1b
				x->parent->color=RED;     // Case 1b
				rightRotate(x->parent, updateNode);  // Case 1b
				w=x->parent->left;       // Case 1b
			}
			if ((w->right->color == BLACK) && (w->left->color == BLACK))
			{
				w->color=RED;  // Case 2b
				x=x->parent;   // Case 2b
			}
			else
			{
				if (w->left->color == BLACK)
				{
					w->right->color=BLACK;  // Case 3b
					w->color=RED;          // Case 3b
					leftRotate(w, updateNode);      // Case 3b
					w=x->parent->left;    // Case 3b
				}

				w->color=x->parent->color;  // Case 4b
				x->parent->color=BLACK;     // Case 4b
				w->left->color=BLACK;      // Case 4b
				rightRotate(x->parent, updateNode);    // Case 4b

				x=root;                  // Case 4b
			}
		}
	} // end of while
	x->color=BLACK;
}

// quite similar to Cormen et al
void RBTree::rbDelete(RBNode *z, void (*updateNode)(RBNode* n, RBTree *T)){
	RBNode *y,*x,*y_oldParent;
	y_oldParent=nil;
	if (z->left == nil || z->right == nil) {
		y = z;  //z has no or one child, deletion is easy
	} else {
		y = treeSuccessor(z);
		y_oldParent=y->parent;
	}

	if (y->left != nil) {
		x=y->left;
	} else {
		x=y->right;
	}
	
	x->parent=y->parent;

	// cut out y:
	if (y->parent == nil) {
		root=x;
	} else {
		if (isLeftChild(y)) {y->parent->left = x;}
			else {y->parent->right = x;}
	}
	
	RBNodecolor old_y = y->color;
	if (y != z) { // 2 children
		//move y to z's old position and delete z
		if (root==z) {
			root=y;
		} else {
			if (isLeftChild(z)) z->parent->left = y;
			else z->parent->right = y;
		}
			
		y->parent = z->parent;
		y->left = z->left;
		y->right = z->right;
		y->color = z->color;  // don't forget to delete z after rbDelete

		y->left->parent = y;
		y->right->parent = y;
		
	}


	if (old_y == BLACK) {
		rbDeleteFixup(x, updateNode);
	}
	
	updateNode(y, this);
	if (y_oldParent!=nil) updateNode(y_oldParent, this);

}


RBNode* RBTree::treeSuccessor(RBNode *x){
	if (x->right != nil) return treeMinimum(x->right);

	RBNode *y = x->parent;
	while ((y!=nil) && (x == y->right)) {
		x=y;
		y=y->parent;
	}
	return y;
}

RBNode* RBTree::treePredecessor(RBNode *x){
	if (x->left != nil) return treeMaximum(x->left);

	RBNode *y = x->parent;
	while ((y!=nil) && (x == y->left)) {
		x=y;
		y=y->parent;
	}
	return y;
}


RBNode* RBTree::treeMinimum(RBNode *x){
	while (x->left != nil) x=x->left;
	return x;
}

RBNode* RBTree::treeMaximum(RBNode *x){
	while (x->right != nil) x=x->right;
	return x;
}


bool RBTree::isLeftChild(RBNode *n){
	#ifndef NDEBUG
	if (n->parent == nil) {
		cerr << "error: isLeftChild, no parent." << endl;
		exit(EXIT_FAILURE);
	}
	#endif
	return (n->parent->left == n);
}

bool RBTree::isRightChild(RBNode *n){
	#ifndef NDEBUG
	if ( n->parent == nil) {
		cerr << "error: isLeftChild, no parent." << endl;
		exit(EXIT_FAILURE);
	}
	#endif
	return (n->parent->right == n);
}


RBNode* RBTree::findRightSiblingLeaf(RBNode *n){
	// go up:
	while (true) {
		if (n->parent!=nil) {
			if (n==n->parent->right)
			  n=n->parent;
			else
			  break;
		} else
			return nil;
	}
		
	n=n->parent;

	// leftmost leaf in n is the right sibling searched
	n=n->right;

	// go down:
	while (n->left != nil) {
		n=n->left;
	}

	return n;
}

RBNode* RBTree::findLeftSiblingLeaf(RBNode *n){
	// go up:
	while (true) {
		if (n->parent!=nil) {
			if (n==n->parent->left)
				n=n->parent;
			else
				break;
		} else
			return nil;
	}
	
	n=n->parent;
	
	// rightmost leaf in n is the left sibling searched
	n=n->left;

	// go down:
	while (n->right != nil) {
		n=n->right;
	}

	return n;
}


int RBTree::getNodeMaxDepth(RBNode *n) {
	int l;
	int r;
	if (n->left == nil) l=0;
	else l=getNodeMaxDepth(n->left);
	if (n->right == nil) r=0;
	else r=getNodeMaxDepth(n->right);

	return (1+(l>r?l:r));
}

int RBTree::getNodeMinDepth(RBNode *n) {
	int l;
	int r;
	if (n->left == 0) l=0;
	else l=getNodeMinDepth(n->left);
	if (n->right == 0) r=0;
	else r=getNodeMinDepth(n->right);

	return (1+(l>r?r:l));
}


void RBTree::printSubTree(RBNode *n){
	cout << ((n->color==BLACK)?"B":"R") << n << " [";
	if (n->left==nil) cout << "N";
		else printSubTree(n->left);
	cout << ",";
	if (n->right==nil) cout << "N";
		else printSubTree(n->right);
	cout << "]";
}


void RBTree::checkSubTree(RBNode *n){

	if (n->left!=nil) {
		if (n->left->parent != n) {
			cout << "au"<< endl;
			exit(1);
		}
		checkSubTree(n->left);
	}

	if (n->right!=nil) {
		if (n->right->parent != n) {
			cout << "au"<< endl;
			exit(1);
		}
		checkSubTree(n->right);
	}

}


void RBTree::checkNode(RBNode *n){
	if (n->left!=nil) {
		if (n->left->parent != n) {
			cout << "au"<< endl;
			exit(1);
		}
		
	}

	if (n->right!=nil) {
		if (n->right->parent != n) {
			cout << "au"<< endl;
			exit(1);
		}
	}	
}

} // namespace

