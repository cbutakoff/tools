#include "MinHeap.h"
#include <iostream>

int main(int argc, char* argv[])
{
    int array[] = {10, 4, 5, 30, 3, 300};

    
    typedef MinHeap<double, double> HeapType;
    std::vector<HeapType::HeapNodeType> vector;
    for(int i=0; i<6; i++)
    {
        HeapType::HeapNodeType node;
        node.key=array[i];
        node.data=array[i]*2;
        vector.push_back(node);
    }
    
    

    HeapType minHeap(vector);

    for(int i=0; i<3; ++i)
    {
        std::cout << minHeap.GetMin() << "  ";
        minHeap.DeleteMin();
    }

    char x;
    std::cin >> x;
 
   return 0;
}
