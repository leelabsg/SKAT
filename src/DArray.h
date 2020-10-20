/*************************************************************
 *
 * NPTest Project
 * File: DArray.h	
 * Date: November 1, 2010
 * Author: Larissa Miropolsky
 *
 * Description:
 *   The Template of Dynamic Array , that includes also 
 *   interator.
 *
 **************************************************************/

#ifndef _DARRAY_H
#define _DARRAY_H



#include <stdio.h>
#include <climits>
#include <cstring>
#include <memory>



/*
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

*/

#pragma once
 
#include <iterator>


#define BYTE char
#define max(a,b)    (((a) > (b)) ? (a) : (b))
#define min(a,b)    (((a) < (b)) ? (a) : (b))



enum CDARRAY_ACTION
{
	NONE,
	REMOVE_RESIZE
};

template<class T> class CDArrayIter;  //ITERATOR


template<class T> class CDArray
{
	
	// Implementation
private:
	
	T** m_pData;		// the actual array of data
	long m_nSize;		// # of elements (upperBound + 1)
	long m_nMaxSize;	// max allocated
	int m_nGrowBy;		// grow amount
public:
	//===========================================================
	// Iterator
	//===========================================================
	friend class CDArrayIter<T>; 
	CDArrayIter<T>* CreateIterator()const{return new CDArrayIter<T>(this);} 
	CDArrayIter<T>* m_iterator;

	T* CurrentItem()
	{
		return m_iterator->currentItem();
	}

	void First()
    {
		this->m_iterator->first();
    }
    void Next()
    {
		this->m_iterator->next();
    }
    bool IsDone()
    {
		return this->m_iterator->isDone();
    }

	//===========================================================
	//  Construction/ Destruction
	//===========================================================
	CDArray()
	{
		m_pData = NULL;
		m_nSize = 0;
		m_nMaxSize = 0;
		m_nGrowBy = 0;
		m_iterator = this->CreateIterator();
	}
	//===========================================================
	CDArray(CDArray* array)
	{
		/*ASSERT(array->GetSize() >= 0);*/
		
		m_pData = NULL;
		m_nSize = 0;
		m_nMaxSize = 0;
		m_nGrowBy = 0;
		
		int inSize = array->GetSize();
		T** inData = array->GetData();
		
		SetSize(inSize);
		memcpy(m_pData, inData, inSize * sizeof(T*));
		m_iterator = this->CreateIterator();
	}
	//===========================================================
	CDArray(const unsigned int size)
	{ 
		m_pData = NULL;
		m_nSize = 0;
		m_nMaxSize = 0;
		m_nGrowBy = 0; 
		SetSize(size); 
		m_iterator = this->CreateIterator();
	}
	//===========================================================
	~CDArray()
	{
		delete [] (BYTE*)m_pData; 
		delete m_iterator;
		m_pData = NULL;
		m_nSize = 0;
		m_nMaxSize = 0;
		m_nGrowBy = 0;
	}
	
	//===========================================================
	// Reset
	//===========================================================
    void Reset()
    {
		// init array slots to NULL
		memset(m_pData, 0, m_nSize*sizeof(T*));
    }
	
	//===========================================================
	// Attributes
	//===========================================================
	long GetSize() { return m_nSize; }   
    bool IsEmpty(){ return GetSize() == 0; }  
	long GetUpperBound() { return m_nSize-1; } 
	
	//===========================================================
	// Clean All
	//===========================================================
	void Free()
	{
		for(long i=0;i<m_nSize;i++)
        {
            //if( IsDuplicated( i ) )
            //    continue;
			if(m_pData[i] != NULL)
			{
				delete m_pData[i];
				m_pData[i] = NULL;
			}
        }
		RemoveAll();
	}
	//===========================================================
	// Clean cell
	//===========================================================
	
	void FreeAt(int nIndex,CDARRAY_ACTION action = NONE)
	{
		if(m_pData[nIndex] != NULL)
		{
			delete m_pData[nIndex];
			m_pData[nIndex] = NULL;
		}
		
		if(action == REMOVE_RESIZE) RemoveAt(nIndex);
	}
	//===========================================================
	// Replace cell
	//===========================================================
	
	void Replace(int nIndex, T* pElement)
	{
		if(m_pData[nIndex] != NULL)
		{
			delete m_pData[nIndex];
			m_pData[nIndex] = NULL;
		}
		SetAt(nIndex,pElement);
	}
	//===========================================================
	// Clean up pointers array
	//===========================================================
	
	void RemoveAll() 
	{ 
		SetSize(0); 
	}
	//===========================================================
	
	void RemoveAt(int nIndex, int nCount = 1)
	{
		/*ASSERT(nIndex >= 0);*/
		/*ASSERT(nCount >= 0);*/
		/*ASSERT(nIndex + nCount <= m_nSize);*/
		
		// just remove a range
		long nMoveCount = m_nSize - (nIndex + nCount);
		
		if (nMoveCount)
			memcpy(&m_pData[nIndex], &m_pData[nIndex + nCount],
			nMoveCount * sizeof(T*));
		m_nSize -= nCount;
	}
	
	//===========================================================
	// Access data
	//===========================================================
	// Accessing elements
	//===========================================================
	T* GetAt(int nIndex)
	{ 
		/*ASSERT(nIndex >= 0 && nIndex < m_nSize);*/
		return m_pData[nIndex]; 
	}
	//===========================================================
    //Last element
	//===========================================================
	T* GetLast(void)
    {
        return GetAt(GetUpperBound());
    }
	//===========================================================
	//First element
	//===========================================================
    T* GetFirst(void)
    {
        return GetAt(0);
    }
	//===========================================================
	// Direct Access to the element data (may return NULL)
	//===========================================================
	T** GetData()	
	{ 
		return m_pData; 
	}
	
	//===========================================================
	// Set data
	//===========================================================
	// SetAt
	//===========================================================
	virtual void SetAt(int nIndex, T* newElement)
	{ 
		/*ASSERT(nIndex >= 0 && nIndex < m_nSize);*/
		m_pData[nIndex] = newElement; 
	}
	//===========================================================
	// Potentially growing the array
	//===========================================================
	virtual void SetAtGrow(int nIndex, T* newElement)
	{
		/*ASSERT(nIndex >= 0);*/
		
		if (nIndex >= m_nSize)
			SetSize(nIndex+1);
		m_pData[nIndex] = newElement;
	}
	//===========================================================
	// Operations that move elements around
	//===========================================================
	virtual void InsertAt(int nIndex, T* newElement, int nCount = 1)
	{
		/*ASSERT(nIndex >= 0);*/    // will expand to meet need
		/*ASSERT(nCount > 0); */    // zero or negative size not allowed
		
		if (nIndex >= m_nSize)
		{
			// adding after the end of the array
			SetSize(nIndex + nCount);  // grow so nIndex is valid
		}
		else
		{
			// inserting in the middle of the array
			long nOldSize = m_nSize;
			SetSize(m_nSize + nCount);  // grow it to new size
			// shift old data up to fill gap
			memmove(&m_pData[nIndex+nCount], &m_pData[nIndex],
				(nOldSize-nIndex) * sizeof(T*));
			
			// re-init slots we copied from
			memset(&m_pData[nIndex], 0, nCount * sizeof(T*));
			
		}
		
		// insert new value in the gap
		/*ASSERT(nIndex + nCount <= m_nSize);*/
		while (nCount--)
			m_pData[nIndex++] = newElement;
	}
	
	//===========================================================
	// Add data
	//===========================================================
	virtual void Init(int val)
	{
		memset(m_pData, val, m_nSize * sizeof(T*));
	}
	//===========================================================
	// Add
	//===========================================================
	virtual int Add(T* newElement)
	{ 
		long nIndex = m_nSize;
		SetAtGrow(nIndex, newElement);
		return nIndex; 
	}
	//===========================================================
	//Append
	//===========================================================
	virtual int Append(CDArray* array)
	{ 
		/*ASSERT(array->GetSize() > 0);*/
		long nIndex = m_nSize;
		
		int inSize = array->GetSize();
		T** inData = array->GetData();
		
		SetSize(m_nSize + inSize);
		memcpy(&m_pData[nIndex], inData,
			inSize * sizeof(T*));
		
		return nIndex; 
	}
	//===========================================================
	//Copy
	//===========================================================
	virtual int Copy(CDArray* array)
	{ 
		/*ASSERT(array->GetSize() > 0);*/
		
		int inSize = array->GetSize();
		T** inData = array->GetData();
		
		SetSize(inSize);
		memcpy(&m_pData[0], inData,
			inSize * sizeof(T*));
		
		return inSize - 1; 
	}
	
	
	//===========================================================
	// Operators
	//===========================================================
	T* operator[](int nIndex) { return GetAt(nIndex); }
	
	
	
	//===========================================================
	// Set size
	//===========================================================
	void SetSize(int nNewSize, int nGrowBy = -1)
	{
		/*ASSERT(nNewSize >= 0);*/
		
		if(nGrowBy != -1)
			m_nGrowBy = nGrowBy;  // set new size
		
		// shrink to nothing
		if(nNewSize == 0)
		{
			delete[] (BYTE*)m_pData;
			m_pData = NULL;
			m_nSize = 0;
			m_nMaxSize = 0;
			m_nGrowBy = 0;
		}
		else 
			// create one with exact size
			if(m_pData == NULL)
			{
				m_pData = (T**) new BYTE[nNewSize * sizeof(T*)];
				memset(m_pData, 0, nNewSize * sizeof(T*));  // zero fill
				m_nSize = nNewSize;
				m_nMaxSize = nNewSize;
			}
			else 
				if(nNewSize <= m_nMaxSize)
				{
					// it fits
					if (nNewSize > m_nSize)
					{
						// initialize the new elements
						memset(&m_pData[m_nSize], 0, (nNewSize-m_nSize) * sizeof(T*));
					}
					m_nSize = nNewSize;
				}
				else
				{
					// otherwise, grow array
					int nGrowBy = m_nGrowBy;
					if (nGrowBy == 0)
					{
						// heuristically determine growth when nGrowBy == 0
						//  (this avoids heap fragmentation in many situations)
						nGrowBy = min(1024, max(4, m_nSize / 8));
					}
					int nNewMax;
					if (nNewSize < m_nMaxSize + nGrowBy)
						nNewMax = m_nMaxSize + nGrowBy;  // granularity
					else
						nNewMax = nNewSize;  // no slush
					
					/*ASSERT(nNewMax >= m_nMaxSize); */ // no wrap around
					
					T** pNewData = (T**) new BYTE[nNewMax * sizeof(T*)];
					
					// copy new data from old
					memcpy(pNewData, m_pData, m_nSize * sizeof(T*));
					
					// construct remaining elements
					/*ASSERT(nNewSize > m_nSize);*/
					
					memset(&pNewData[m_nSize], 0, (nNewSize-m_nSize) * sizeof(T*));
					
					// get rid of old stuff (note: no destructors called)
					delete[] (BYTE*)m_pData;
					m_pData = pNewData;
					m_nSize = nNewSize;
					m_nMaxSize = nNewMax;
				}
	}

};



template<class T> class CDArrayIter 
{
private:
	friend class CDArray<T>;
	const CDArray<T>* cdarr;
	int index;

public:
    CDArrayIter(const CDArray<T> *s)
    {
        cdarr = s;
		index = 0;
    }
    void first()
    {
        index = 0;
    }
    void next()
    {
        index++;
    }
    bool isDone()
    {
        return index == cdarr->m_nSize;
    }
    T* currentItem()
    {
        return cdarr->m_pData[index];
    }
};

#endif // _DARRAY_HPP
