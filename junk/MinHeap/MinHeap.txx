/*=========================================================================
Copyright (c)  Anton Kaminsky, http://www.codeproject.com/Tips/816934/Min-Binary-Heap-Implementation-in-Cplusplus
Adapted by Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
#ifndef MinHeap_txx
#define MinHeap_txx

#include "MinHeap.h"
#include <iostream>

template<class StorageDataType, class KeyType>
MinHeap<StorageDataType,KeyType>::MinHeap(const std::vector<HeapNodeType>& vector) : _vector(vector)
{
    Heapify();
}

template<class StorageDataType, class KeyType>
MinHeap<StorageDataType,KeyType>::MinHeap() 
{
}

template<class StorageDataType, class KeyType>
void MinHeap<StorageDataType,KeyType>::Heapify()
{
    int length = _vector.size();
    for(int i=length-1; i>=0; --i)
    {
        BubbleDown(i);
    }
}

template<class StorageDataType, class KeyType>
void MinHeap<StorageDataType,KeyType>::BubbleDown(int index)
{
    int length = _vector.size();
    int leftChildIndex = 2*index + 1;
    int rightChildIndex = 2*index + 2;

    if(leftChildIndex >= length)
        return; //index is a leaf

    int minIndex = index;

    if(_vector[index].key > _vector[leftChildIndex].key)
    {
        minIndex = leftChildIndex;
    }
    
    if((rightChildIndex < length) && (_vector[minIndex].key > _vector[rightChildIndex].key))
    {
        minIndex = rightChildIndex;
    }

    if(minIndex != index)
    {
        //need to swap
        HeapNodeType temp = _vector[index];
        _vector[index] = _vector[minIndex];
        _vector[minIndex] = temp;
        BubbleDown(minIndex);
    }
}

template<class StorageDataType, class KeyType>
void MinHeap<StorageDataType,KeyType>::BubbleUp(int index)
{
    if(index == 0)
        return;

    int parentIndex = (index-1)/2;

    if(_vector[parentIndex].key > _vector[index].key)
    {
        HeapNodeType temp = _vector[parentIndex];
        _vector[parentIndex] = _vector[index];
        _vector[index] = temp;
        BubbleUp(parentIndex);
    }
}

template<class StorageDataType, class KeyType>
void MinHeap<StorageDataType,KeyType>::Insert(KeyType key, StorageDataType& data)
{
    HeapNodeType node;
    node.key = key;
    node.data = data;

    
    int length = _vector.size();
    //_vector[length] =  node;
    _vector.push_back(node);
    
    BubbleUp(length);
}

template<class StorageDataType, class KeyType>
StorageDataType MinHeap<StorageDataType,KeyType>::GetMin()
{
    return _vector[0].data;
}
    
template<class StorageDataType, class KeyType>
void MinHeap<StorageDataType,KeyType>::DeleteMin()
{
    int length = _vector.size();

    if(length == 0)
    {
        return;
    }
    
    _vector[0] = _vector[length-1];
    _vector.pop_back();

    BubbleDown(0);
}

#endif
