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

// ------ Dynamic Sequence with Indels -----
// uses huffman shaped wavelet tree
// space: O(n H_0) time: O((log n) H_0)
// papers: V. Maekinen, G. Navarro. Dynamic Entropy-Compressed Sequences and Full-Text
//           Indexes. CPM 2006, Chapter 3.6 

#ifndef GUARD_DynFMI
#define GUARD_DynFMI


#include <iostream>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <bitset>
#include <stdio.h>
#include <string.h>

#ifndef SAMPLE
#define SAMPLE 1024
#endif

#ifndef BUFFER
#define BUFFER 1048576
#endif


#include "bvtree.h"
#include "handle.h"
#include "pos.h"

#if SAMPLE!=0
#include "ssatree.h"
#endif


#ifndef uchar
#define uchar unsigned char
#endif


namespace dynfmi { class DynFMI; }

using dynfmi::DynFMI;

namespace dynfmi {

using bvtree::BVTree;
using postree::PosTree;
using handletree::HandleTree;
using sampledsatree::SampledSATree;


//TODO in Klasse schieben
const size_t sampleInterval = SAMPLE;


class WaveletNode{
	public:
	
	WaveletNode *left;
	WaveletNode *right;
	WaveletNode *parent;
	float relativeFrequency; // used only while construction
	uchar c0;      // used also while construction
	uchar c1;
	
	BVTree *bittree;
	
	WaveletNode(uchar c, float relativeFrequency)
		: left(0), right(0), parent(0), relativeFrequency(relativeFrequency), c0(c), bittree(0){}

	WaveletNode(WaveletNode *left, WaveletNode *right)
		: left(left), right(right), parent(0), bittree(0) {
		relativeFrequency = left->relativeFrequency + right->relativeFrequency;
		left->parent = this;
		right->parent = this;
	}
	
	~WaveletNode(){
		delete bittree;
	}
	
	bool operator>(const WaveletNode &a) const {
		return (relativeFrequency > a.relativeFrequency);
	}
	
	
};

  
class DynFMI{
	public:
		
				
		void initEmptyDynFMI(const float *f); //argument: array length 256 containing relative frequencies of characters
		void empty();
		
		
		
		
		/**
		You can use these functions to compute the relative character frequencies.
		Example:  createUniformCharDistribution("ACGT", 10000000, 500) - means 500 DNA-sequences
		of overall length  10000000; including separation character, probabilities are 1/5.
		**/
		static float* createUniformCharDistribution(const uchar *sampleText, size_t expectedTotalLength, size_t expectedNumberOfTexts); //for equal length encoding of characters
		static float* createCharDistribution(const uchar *sampleText, size_t expectedTotalLength, size_t expectedNumberOfTexts); //for variable length encoding of characters
		
		
		size_t addText(const uchar *str, size_t n); // length n of text must include "\0"!
		size_t addTextFromFile(char* filename, size_t n); //data is read directly from disk using only BUFFER(macro) bytes
		uchar* retrieveText(size_t text);
		bool deleteText(size_t key);
		
		
		size_t count(const uchar *pattern); // only counting, no locate!
		size_t* locate(const uchar *pattern); // needs macro SAMPLE > 0! (default)
		// returns array: [length,{handle, position,}^length]
		
		pair<size_t, size_t>* backwardssearch(const uchar *pattern);
		// returns pair: (lower bound, upper bound)
		pair<size_t, size_t>* SAlookup(size_t j); // needs macro SAMPLE > 0! (default)
		// returns pair: (handle, matchposition)

		//LF(i)-mapping: C[L[i]]+rank_L[i](L,i)
		size_t LFmapping(size_t i) {
			uchar s=(*this)[i];
			return getNumberOfSymbolsSmallerThan(s) + rank(s,i);
			}
		
		
		size_t getSize() {return root->bittree->getPositions();}
		size_t getCollectionSize() { return pos->getSize();}
		
		uchar*   getBWT();
		ostream& getBWTStream(ostream& stream);

		uchar operator[](size_t i);
		size_t* getKeys() {return handle->getKeys(); }	
		

		
		~DynFMI();

	private:
		WaveletNode *root;    // root of the wavelet tree
		WaveletNode **leaves; // needed for construction and select
		
		size_t codes[256];
		int codelengths[256];
		size_t C[256+256];
		
#if SAMPLE!=0
		BVTree *SampledSAPositionsIndicator;
#endif		
		size_t iterate;
				
		HandleTree *handle;
		PosTree *pos;
#if SAMPLE!=0
		SampledSATree *sampledSATree;
#endif


		size_t rank(uchar c, size_t i);
		size_t select(uchar c, size_t i);
	
		void insert(uchar c, size_t i);
		void deleteSymbol(size_t i);



		// functions
		size_t getNumberOfSymbolsSmallerThan(uchar c);
		
		
		
		void makeCodes(size_t code, int bits, WaveletNode *node);
		void deleteLeaves(WaveletNode *node);
		void appendBVTrees(WaveletNode *node);
		
		void deleteDynFMINodes(WaveletNode *n);
		
		//Iterator
		void iterateReset();
		bool iterateNext();
		uchar iterateGetSymbol();
		void recursiveIterateResetWaveletNode(WaveletNode *w);
		

		// small help functions
		static double log2(double x){
			return (log10(x) / log10((double)2));
		}
		
		static int binaryTree_parent(int i){
			return (int)floor((double)i/2);
		}

		static int binaryTree_left(int i){
			return 2*i;
		}

		static int binaryTree_right(int i){
			return 2*i + 1;
		}

		static bool binaryTree_isLeftChild(int i){
			return (i%2==(int)0);
		}

		static bool binaryTree_isRightChild(int i){
			return (i%2==(int)1);
		}
		
};

} //namespace



namespace std
{

template<> struct greater<dynfmi::WaveletNode*>
  {
    bool operator()(dynfmi::WaveletNode const* p1, dynfmi::WaveletNode const* p2)
    {
      if(!p1)
        return false;
      if(!p2)
        return true;
      return *p1 > *p2;
    }
  };
}


#endif

