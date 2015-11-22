/*=========================================================================
Copyright (c)  Anton Kaminsky, http://www.codeproject.com/Tips/816934/Min-Binary-Heap-Implementation-in-Cplusplus
Adapted by Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
#ifndef MinHeap_h
#define MinHeap_h

#include <vector>
#include <memory>


template<class StorageDataType, class KeyType>
class MinHeap
{
public:
    typedef struct __fnode
    {
        KeyType key; //sorting key
        StorageDataType data; //data to be stored
    } HeapNodeType;


    MinHeap(const std::vector<HeapNodeType>& vector);
    MinHeap();

    void Insert(KeyType key, StorageDataType& data);
    StorageDataType GetMin();
    void DeleteMin();
    
private:
    std::vector<HeapNodeType> _vector;
    void BubbleDown(int index);
    void BubbleUp(int index);
    void Heapify();
    
};


#include "MinHeap.txx"

#endif
