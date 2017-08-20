//
// Created by tyz on 2015/11/24.
//

#ifndef PPR_HEAP_H
#define PPR_HEAP_H


/*
 * C++ Program to Implement Binary Heap
 */
#include <iostream>
#include <cstdlib>
#include <vector>
#include <iterator>

using namespace std;

/*
 * Class Declaration
 */
template<typename T, typename Compare = less<T>>
class BinaryHeap {
private:
    vector<pair<T, int>> heap;
    vector<int> heap_idx;

    int left(int parent) {
        int l = 2 * parent + 1;
        if (l < int(heap.size()))
            return l;
        else
            return -1;
    }


    int right(int parent) {
        int r = 2 * parent + 2;
        if (r < int(heap.size()))
            return r;
        else
            return -1;
    }

    int parent(int child) {
        int p = (child - 1) / 2;
        if (child == 0)
            return -1;
        else
            return p;
    }


    void heapswap(int x, int y) {
        swap(heap_idx[heap[x].second], heap_idx[heap[y].second]);
        swap(heap[x], heap[y]);
    }

    void heapifyup(int in) {
        while (true) {
            int pin = parent(in);
            // if smaller than parent then swap
            if (in >= 0 && pin >= 0 && comp(heap[in].first, heap[pin].first)) {
                heapswap(in, pin);
                in = pin;
            }
            else
                break;
        }
    }


    void heapifydown(int in) {
        while (true) {
            int child = left(in);
            int child1 = right(in);
            // child point to smaller one
            if (child >= 0 && child1 >= 0 && comp(heap[child1].first, heap[child].first)) {
                child = child1;
            }
            // if child is smaller, then swap
            if (child > 0 && comp(heap[child].first, heap[in].first)) {
                heapswap(in, child);
                in = child;
            }
            else {
                break;
            }
        }
    }

    Compare comp;
    // comp return true means more close to root
public:
    BinaryHeap(int maxn, Compare _comp = std::less<T>()) {
        heap_idx = vector<int>((unsigned long) maxn + 10, -1);
        comp = _comp;
//        cout << comp(1, 2) << endl;
//        cout << comp(2, 1) << endl;
    }

    map<int, T> as_map() {
        map<int, T> rtn;
        for (auto item: heap) {
            rtn[item.second] = item.first;
        }
        return rtn;
    };

    vector<pair<T, int>>& get_elements(){
        return heap;
    }

    T get_value(int idx) {
        return heap[heap_idx[idx]].first;
    };

    void insert(int index, T element) {
        heap_idx[index] = (int) heap.size();
        heap.push_back(MP(element, index));
        heapifyup(int(heap.size()) - 1);
    }

    void delete_top() {
        if (heap.size() == 0) {
            cout << "Heap is Empty" << endl;
            return;
        }
        heapswap(0, (int) (heap.size() - 1));
        heap_idx[heap.back().second] = -1;
        heap.pop_back();
        heapifydown(0);
        //cout << "Element Deleted" << endl;
    }

    void modify(int index, T newvalue) {
        int idx = heap_idx[index];
        if (comp(newvalue, heap[idx].first)) {
            heap[idx].first = newvalue;
            heapifyup(idx);
        }
        else {
            heap[idx].first = newvalue;
            heapifydown(idx);
        }
    }

    pair<T, int> extract_top() {
        return heap.front();
    }

    void display() {
        cout << "Heap -->  ";
        for (auto p:heap) {
            cout << p.first << " ";
        }
        cout << endl;

    }

    int size() {
        return int(heap.size());
    }

    void clear() {
//        heap.clear();
        while (size())
            delete_top();
    }

    // return true is has index
    bool has_idx(int idx) {
        return heap_idx[idx] != -1;
    }

    void verify() {
        for (int child = 0; child < heap.size(); ++child) {
            int p = parent(child);
            if (p >= 0)
                // child cannot be smaller than parent
                assert (!comp(heap[child].first, heap[p].first));

        }
    }
};


#endif //PPR_HEAP_H
