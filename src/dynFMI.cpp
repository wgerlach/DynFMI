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
#include <fstream>
#include <stack>
#include <queue>
#include <functional>
#include <algorithm>

#include "dynFMI.h"






// -------- DynFMI --------
namespace dynfmi {

DynFMI::~DynFMI(){
	empty();
}

void DynFMI::empty(){
	deleteDynFMINodes(root);
	delete[] leaves;
#if SAMPLE!=0
	delete SampledSAPositionsIndicator;
#endif
	delete handle;
	delete pos;
#if SAMPLE!=0
	delete sampledSATree;
#endif
}


void DynFMI::deleteDynFMINodes(WaveletNode *n){
	if (n->right) deleteDynFMINodes(n->right);
	if (n->left) deleteDynFMINodes(n->left);
	
	delete n;
}


void DynFMI::iterateReset(){
	iterate = 1;
	recursiveIterateResetWaveletNode(root);
}

void DynFMI::recursiveIterateResetWaveletNode(WaveletNode *w){
	w->bittree->iterateReset();
	
	if (w->left) recursiveIterateResetWaveletNode(w->left);
	if (w->right) recursiveIterateResetWaveletNode(w->right);
}


bool DynFMI::iterateNext(){
	iterate++;
	return !(iterate > getSize());
}

uchar DynFMI::iterateGetSymbol(){

	size_t i = iterate;
	WaveletNode *walk = root;	
	bool bit;
		
	while (true) {
		
		bit = walk->bittree->iterateGetBit(); // TODO improve
		i=walk->bittree->iterateGetRank(bit);
		
		walk->bittree->iterateNext();
		
		
		if (bit) { //bit = 1
			if (walk->right == 0) return walk->c1;
			walk=walk->right;
		} else { // bit = 0
			if (walk->left == 0) return walk->c0;
			walk=walk->left;
		}
		
		
	} // end of while
	
}


uchar* DynFMI::getBWT(){
	size_t length = root->bittree->getPositions();
	length++;
	uchar *text = new uchar[length];
	bool data=true;
	
	size_t i = 0;
	
	iterateReset();
	
	while (data) {
		text[i] = iterateGetSymbol();	
		data = iterateNext();
		i++;
	}
	text[length-1]=0;
	
	return text;
}

ostream& DynFMI::getBWTStream(ostream& stream){
	bool data=true;
	
	iterateReset();
	
	while (data) {
		stream << iterateGetSymbol();
		data = iterateNext();
	}

	return stream;
}

//void DynFMI::printDynFMIContent(ostream& stream){
//	uchar c;
//	for (size_t i=1; i<=getSize(); i++) 
//	{
//		c =(*this)[i];
//		if (c==0) c= '#';
//		stream << c;
//	}


void DynFMI::deleteLeaves(WaveletNode *node){
	bool leaf = true;

	if (node->left) {
		// internal node
		leaf = false;
		deleteLeaves(node->left);
		
	}
	if (node->right){
		leaf = false;
		deleteLeaves(node->right);
	} 
	
	if (leaf) {
		// is a leaf, delete it!
		if (node->parent) {
			if (node==node->parent->left) node->parent->left=0;
				else node->parent->right=0;
		}
		delete node;
	}
}

void DynFMI::makeCodes(size_t code, int bits, WaveletNode *node){
	#ifndef NDEBUG
	if (node == node->left) {
		cout << "makeCodes: error" << endl;
		exit(0);
		}
	#endif

	if (node->left) {
		makeCodes(code | (0 << bits), bits+1, node->left);
		makeCodes(code | (1 << bits), bits+1, node->right);
	} else {
		codes[node->c0] = code;
		codelengths[node->c0] = bits+1;
	}
}

void DynFMI::appendBVTrees(WaveletNode *node){
	node->bittree = new BVTree();

	if (node->left) appendBVTrees(node->left);
	if (node->right) appendBVTrees(node->right);
}


void DynFMI::initEmptyDynFMI(const float *f){
	// pointers to the leaves for select
	leaves = (WaveletNode**) new WaveletNode*[256];
	for(int j=0; j<256; j++) leaves[j]=0; 


	
#ifndef NDEBUG
	float sum = 0;
	for (int i=0; i<256; i++) {
		if (f[i] < 0) {
			cerr << "initEmptyDynFMI: f[i] must not be negative! i=" << i << endl;
			exit(0);
		}
		sum += f[i];
	}
	if ((sum < 0.95) || (sum > 1.05)) {
		cerr << "sum != 1: " << sum << endl;
		exit(0);
	} else  cout << "initEmptyDynFMI, ok: sum= " << sum << endl;
#endif

	for(int j=0; j<256; j++) { //all possible characters
		if (f[j]!=0) { //only those that exist
			leaves[j] = new WaveletNode((uchar)j, f[j]); 
		}
	}
	
	
	// Huffman shape, Veli's approach:
	priority_queue< WaveletNode*, vector<WaveletNode*>, greater<WaveletNode*> > q;
	
	
	for(int j=0; j<256; j++){
		if (leaves[j]!=0) {
			q.push( (leaves[j]) );
		}
		codes[j] = 0;
		codelengths[j] = 0;
	}
	
	// creates huffman shape:
	while (q.size() > 1) {
		
		WaveletNode *left = q.top();
		q.pop();
		
		WaveletNode *right = q.top();
		q.pop();
		
		q.push(  new WaveletNode(left, right) );
	}	
	

	root = q.top();
	q.pop();
	
			
	makeCodes(0,0, root);	// writes codes and codelengths

#ifndef NDEBUG
	for(int j=0; j<256; j++){
	cout << "code " << j << ": " << codes[j]<< endl;
	}
#endif
	
	// merge leaves	(one leaf represent two characters!)
	for(int j=0; j<256; j++){
	
		if (leaves[j]) {
		
			if (leaves[j]->parent->left==leaves[j]) {
				leaves[j]->parent->c0=j;
			} else {
				leaves[j]->parent->c1=j;
			}
			leaves[j]=leaves[j]->parent; // merge
		}
	}

	
	deleteLeaves(root);
	
	appendBVTrees(root);
	
	// array C needed for backwards search
	for(int j=0; j<256+256; j++) C[j] = 0;
	
#if SAMPLE!=0
	this->SampledSAPositionsIndicator = new BVTree();
	this->sampledSATree = new SampledSATree();
#endif
	
	this->pos = new PosTree(sampleInterval);
	this->handle = new HandleTree();
	
	pos->handle= this->handle;
	handle->pos=this->pos;
	
}



void DynFMI::insert(uchar c, size_t i){
	#ifndef NDEBUG
	if (codelengths[c]==0) {
		cerr << "error: Symbol \"" << c << "\" (" << (int)c << ") is not in the code table!" << endl;
		exit(EXIT_FAILURE);
	}
	#endif
	
	size_t level = 0;
	size_t code = codes[c];

	bool bit;
	WaveletNode *walk = root;	
		
	while (walk) {
		
		bit = ((code & (1u << level)) != 0); 
		
		walk->bittree->insertBit(bit,i); // TODO improve
		i=walk->bittree->rank(bit, i);

		if (bit) { //bit = 1
			walk=walk->right;
		} else { // bit = 0
			walk=walk->left;
		}
		
		level++;		
	} // end of while
	
	int j = 256+c;
	while(j>1) {
		C[j]++;
		j=binaryTree_parent(j);
		}
	C[j]++;	
	
}

void DynFMI::deleteSymbol(size_t i){
	WaveletNode *walk = root;	
	bool bit;
	uchar c;
	while (true) {

		// original slow version:
		//bit = (*walk->bittree)[i];
		//old_i = i;		
		//i=walk->bittree->rank(bit, i);	
		//walk->bittree->deleteBit(old_i);


		walk->bittree->deleteBit(i);
		bit=walk->bittree->getLastBitDeleted();
		i=walk->bittree->getLastRank();
		
		if (bit) { //bit = 1
			if (walk->right == 0) {
				c = walk->c1;
				break;
				}
			walk=walk->right;
		} else { // bit = 0
			if (walk->left == 0) {
				c = walk->c0;
				break;
			}
			walk=walk->left;
		}
		
		
	} // end of while
	
	int j = 256+c;
	while(j>1) {
		C[j]--;
		j=binaryTree_parent(j);
		}
	C[j]--;	
	
	
}


uchar DynFMI::operator[](size_t i){
	WaveletNode *walk = root;	
	bool bit;
		
	while (true) {
		
		bit = (*walk->bittree)[i]; //TODO improve by reducing
		i=walk->bittree->rank(bit, i);

		if (bit) { //bit = 1
			if (walk->right == 0) return walk->c1;
			walk=walk->right;
		} else { // bit = 0
			if (walk->left == 0) return walk->c0;
			walk=walk->left;
		}
		
		
	} // end of while
}

size_t DynFMI::rank(uchar c, size_t i){
	if (i==0) return 0;

	size_t level = 0;
	size_t code = codes[c];
	
	
	bool bit;
	WaveletNode *walk = root;	
		
	while (true) {
		
		bit = ((code & (1u << level)) != 0);
		
		i=walk->bittree->rank(bit, i);
		if (bit) { //bit = 1
			if (walk->right == 0) return i;
			walk=walk->right;
		} else { // bit = 0
			if (walk->left == 0) return i;
			walk=walk->left;
		}
	
		level++;		
	} // end of while
	
}

size_t DynFMI::select(uchar c, size_t i){
	
	WaveletNode *walk = leaves[c];	
	
	bool bit = (walk->c1==c);
	
	while (walk->parent) {
		i=walk->bittree->select(bit, i);
		
		bit = (walk == walk->parent->right);
		walk=walk->parent;
	} // end of while
	
	i=walk->bittree->select(bit, i);

	return i;
}

pair<size_t, size_t>* DynFMI::backwardssearch(const uchar *pattern){

	size_t i = strlen((char*)pattern);
	size_t sp = 1;

	size_t ep = getSize();
	uchar s;
	while ((sp <= ep) && (i >=1)) {
		s=pattern[i-1];
	//cout << s << "" << getNumberOfSymbolsSmallerThan(s) << " and " << rank(s, sp-1UL) << endl;
		sp=getNumberOfSymbolsSmallerThan(s) + rank(s, sp-1UL) + 1UL;
		ep=getNumberOfSymbolsSmallerThan(s) + rank(s,ep);
		
		i--;
	}

	pair<size_t, size_t> *p = new pair<size_t, size_t>(sp, ep);
	return p;
}



size_t DynFMI::count(const uchar *pattern){

	pair<size_t, size_t> *p = backwardssearch(pattern);

	return p->second - p->first +1;
}

#if SAMPLE!=0
size_t* DynFMI::locate(const uchar *pattern){

	size_t i;
	size_t sp, ep;
	

	pair<size_t, size_t> *p = backwardssearch(pattern);
	sp = p->first;
	ep = p->second;


	size_t numberOfMatches = ep-sp +1;
	size_t *matches = new size_t[2*numberOfMatches + 1];
	matches[0] = numberOfMatches;

	size_t diff;
	size_t k = 1;
	uchar s;

	for (i=sp; i <= ep; i++) { // each match

		diff = 0;
		while (!(*SampledSAPositionsIndicator)[i]) {
			diff++;
			s=(*this)[i];
			i=getNumberOfSymbolsSmallerThan(s) + rank(s, i);
	
		}
	

		sampledSATree->getSample(SampledSAPositionsIndicator->rank(true, i));

		//cout << "found: (" << sampledSATree->handle  << ", " << sampledSATree->sampledSAValue + diff << ")"<< endl;
		matches[k++] = sampledSATree->handle;
		matches[k++] = sampledSATree->sampledSAValue + diff;
	}

	return matches;
}
#else
size_t* DynFMI::locate(const uchar *pattern){ // if SAMPLE ist turned off
	return 0;
}
#endif

#if SAMPLE!=0
pair<size_t, size_t>* DynFMI::SAlookup(size_t i){
		size_t diff = 0;
		uchar s;
		while (!(*SampledSAPositionsIndicator)[i]) {
			diff++;
			s=(*this)[i];
			i=getNumberOfSymbolsSmallerThan(s) + rank(s, i);
	
		}
	

		sampledSATree->getSample(SampledSAPositionsIndicator->rank(true, i));

		//cout << "found: (" << sampledSATree->handle  << ", " << sampledSATree->sampledSAValue + diff << ")"<< endl;

		//Match *match = new Match (sampledSATree->handle, sampledSATree->sampledSAValue + diff);
		pair<size_t, size_t> *match = new pair<size_t, size_t>(sampledSATree->handle, sampledSATree->sampledSAValue + diff);

		return match;
}
#else
pair<size_t, size_t>* DynFMI::getSAValue(size_t j){ // if SAMPLE ist turned off
		
		return 0;		
}
#endif

// length n of text must include "\0"!
size_t DynFMI::addText(const uchar *str, size_t n){
	if (!root) return 0;

	size_t i;
#if SAMPLE!=0
	bool sample = false;
#endif
	i=pos->getSize()+1;
	
	size_t key = pos->appendText(n);

	insert(str[n-2],i); // insert second last character, corresponds to suffix of length 1
#if SAMPLE!=0
	sample = (n-1 % sampleInterval == 0);
	SampledSAPositionsIndicator->insertBit(sample, i);	
	if (sample) sampledSATree->insertSample(n-1, key, SampledSAPositionsIndicator->rank(true,i));
#endif
	for (size_t t=n-2; t > 0; t--) {
		i=1UL+getNumberOfSymbolsSmallerThan(str[t]) + rank(str[t],i);
		insert(str[t-1],i);
#if SAMPLE!=0
		sample = ((t) % sampleInterval == 0);
		SampledSAPositionsIndicator->insertBit(sample, i);

		if (sample) sampledSATree->insertSample(t, key, SampledSAPositionsIndicator->rank(true,i));
#endif
	}
	

	i=1UL + getNumberOfSymbolsSmallerThan(str[0]) + rank(str[0],i);
	insert(str[n-1],i);
#if SAMPLE!=0
	SampledSAPositionsIndicator->insertBit(true, i); 
	sampledSATree->insertSample(1, key, SampledSAPositionsIndicator->rank(true,i));
#endif	
	
	return key;
}


size_t DynFMI::addTextFromFile(char* filename, size_t n){
	size_t i;
#if SAMPLE!=0
	bool sample = false;
#endif
	i=pos->getSize()+1;


	std::ifstream str(filename, ios::binary);
	if (!str)
	{
		cerr << "error reading file: " << filename << endl;
		exit(EXIT_FAILURE);
	}

	char c;
	//char oldc;


	size_t buf_len; //1048576; //10MB
	char * buffer;
	size_t buf_pos;

	if (BUFFER==0) buf_len = n;
		else buf_len = BUFFER;

	buffer = new char[buf_len];
	buf_pos = (n-2) / buf_len;
	str.seekg(buf_pos*buf_len);
	str.read(buffer, buf_len);
	str.clear(std::ios_base::goodbit);


	size_t key = pos->appendText(n);
	
	//str.seekg(n-2);

	//c=str.get();
	c=buffer[(n-2)%buf_len];
	insert(c,i); // insert second last character, corresponds to suffix of length 1
#if SAMPLE!=0
	sample = (n-1 % sampleInterval == 0);
	SampledSAPositionsIndicator->insertBit(sample, i);	
	if (sample) sampledSATree->insertSample(n-1, key, SampledSAPositionsIndicator->rank(true,i));
#endif
	if ((n-2)%buf_len == 0) {
		c=buffer[0];
		buf_pos--;
		str.seekg(buf_pos*buf_len);
		str.read(buffer, buf_len);
	}
	for (size_t t=n-2; t > 0; t--) {

		i= 1UL + getNumberOfSymbolsSmallerThan(c) + rank(c,i);

		c=buffer[(t-1)%buf_len];
		insert(c,i);

		//check buffer
		if (((t-1)%buf_len == 0) && ((t-1)!=0)) {
#ifndef NDEBUG
			if (buf_pos == 0) {
				cerr << "buf_pos too small" << endl;
				exit(0);
			}
#endif
			buf_pos--;
			str.seekg(buf_pos*buf_len);
			str.read(buffer, buf_len);
		}

#if SAMPLE!=0
		sample = ((t) % sampleInterval == 0);
		SampledSAPositionsIndicator->insertBit(sample, i);

		if (sample) sampledSATree->insertSample(t, key, SampledSAPositionsIndicator->rank(true,i));
#endif
	}

	i= 1UL + getNumberOfSymbolsSmallerThan(c) + rank(c,i);
	insert(c,i);

#if SAMPLE!=0
	SampledSAPositionsIndicator->insertBit(true, i); 
	sampledSATree->insertSample(1, key, SampledSAPositionsIndicator->rank(true,i));
#endif	

	str.close();

	return key;
}


uchar* DynFMI::retrieveText(size_t key){
	#ifndef NDEBUG		
	cout << "key: " << key << endl;
	#endif
	uchar *text=0;
	// TODO access two times, too much
	size_t i=handle->getPos(key); //=bwtEnd;
	if (i == 0 ) {
		return text;
	}

	size_t n=pos->getTextSize(i);
	
	text = new uchar[n]; // last byte 0 for cout
	
	uchar c;
	
//TODO better	
	for (size_t t=n-2; t < n; t--) {
		c = (*this)[i];
		text[t]=(c==0?'#':c);
		
		i= getNumberOfSymbolsSmallerThan(c) + rank(c,i);
	}

	#ifndef NDEBUG		
	for (size_t i=0; i< n-1 ; i++) if (text[i]==0) text[i]='-';
	#endif	
	
	text[n-1] = 0;
	return text;
}


bool DynFMI::deleteText(size_t key){
	// TODO access two times, can do better
	size_t i=handle->getPos(key); //=bwtEnd;
	if (i == 0) return false;

	size_t n=pos->getTextSize(i);
	size_t *text = new size_t[n]; 
	uchar c;
	
	for (size_t t=n-1; t <n; t--) {
		c = (*this)[i];
		text[t]=i;
		i= getNumberOfSymbolsSmallerThan(c) + rank(c,i); // TODO improve by lastrank!
	}


	sort(text, text + n);


	for (i=n-1; i<n; i--) {
		this->deleteSymbol(text[i]); 
#if SAMPLE!=0
		SampledSAPositionsIndicator->deleteBit(text[i]);		
//TODO samples ?
#endif
	}
	delete[] text;
	handle->deleteKey(key);
	
		
	return true;
}

size_t DynFMI::getNumberOfSymbolsSmallerThan(uchar c){
	int j = 256+c;
	size_t r=0;
	while(j>1) {
		if (binaryTree_isRightChild(j)) 
			r += C[binaryTree_left(binaryTree_parent(j))];
		
		j=binaryTree_parent(j);
	}
	return r;
}






float* DynFMI::createUniformCharDistribution(const uchar *sampleText, size_t expectedTotalLength, size_t expectedNumberOfTexts){
	bool count[256];
	count['\0'] = true;
	for (int i=1; i<256; i++) count[i] = false;

	size_t i=0;
	size_t sampleLength;

	while (sampleText[i] != '\0') {
		count[sampleText[i]] = true;
		i++;
	}

	sampleLength=i;

	

	int alphabetSize = 0;
	
	for (int i=0; i<256; i++) if (count[i]) alphabetSize++;
	

	float *f = new float[256];
	for (int i=0; i<256; i++) f[i]=0;
	

	for (int i=0; i<256; i++) {
		if (count[i]) f[i] = 1/(double)alphabetSize; 
#ifndef NDEBUG	
		if ((i >= 39) && (i<=122)) cout << (char)i << ": " << f[i] << endl;
		else cout << i << ": " << f[i] << endl;
#endif		

	}


#ifndef NDEBUG
	float sum = 0;
	for (int i=0; i<256; i++) {
		if (f[i] < 0) {
			cerr << "createUniformCharDistribution: f[i] must not be negative! i=" << i << endl;
			exit(0);
		}

		sum += f[i];
	}
	if ((sum < 0.95) || (sum > 1.05)) {
		cerr << "createUniformCharDistribution: sum != 1: " << sum << endl;
		exit(0);
	} else  cout << "createUniformCharDistribution, ok: uniform sum= " << sum << endl;
#endif


	return f;
}

float* DynFMI::createCharDistribution(const uchar *sampleText, size_t expectedTotalLength, size_t expectedNumberOfTexts){
	size_t count[256];

	for (int i=0; i<256; i++) count[i] = 0;

	size_t i=0;
	size_t sampleLength;

	while (sampleText[i] != '\0') {
		count[sampleText[i]]++;
		i++;
	}

	//for (int i=0; i<256; i++) cout << count[i] << endl;

	sampleLength=i;

	for (int i=1; i<256; i++) {
		if (count[i] != 0) count[i] = (size_t)(count[i]*(expectedTotalLength/(double)sampleLength));  // scale
	}
	count['\0'] = expectedNumberOfTexts;

	


#ifndef NDEBUG
	for (int i=0; i<256; i++) cout << count[i] << endl;

	if (count['\0'] == 0) {
		cerr << "count['\0'] = 0" << endl;
		exit(0);
	}
#endif

	float *f = new float[256];

	//sampleLength += count['\0'];

	for (int i=0; i<256; i++) {
		if (count[i]) {
			f[i] = (float)(count[i]/(double)(expectedTotalLength + expectedNumberOfTexts));
			}
		else {
			f[i]=0;
		}
#ifndef NDEBUG	
	if ((i >= 39) && (i<=122)) cout << (char)i << ": " << f[i] << endl;
	else cout << i << ": " << f[i] << endl;
#endif		
	}


#ifndef NDEBUG
	float sum = 0;
	for (int i=0; i<256; i++) {
		if (f[i] < 0) {
			cerr << "createCharDistribution: f[i] must not be negative! i=" << i << endl;
			exit(0);
		}

		sum += f[i];
	}
	if ((sum < 0.95) || (sum > 1.05)) {
		cerr << "createCharDistribution: sum != 1: " << sum << endl;
		exit(0);
	} else  cout << "createCharDistribution, ok: sum= " << sum << endl;
#endif

	return f;
}


} // namespace

