/*
 * File:   Node.hpp
 * Author: Yuri Meiburg <yuri@meiburg.nl>
 *
 * Created on July 7, 2011, 6:00 PM
 */

#ifndef NODE_HPP
#define NODE_HPP
using namespace std;
template< class T > class Node {
  private:
    T value;
    vector< Node< T > * > children;
  public:
    Node(T v): value(v) {}; //{ value = v; } ;
    ~Node() {
        typename vector< Node< T >* >::iterator n = children.begin();
        while (n != children.end()) {
            children.erase(n);
            n = children.begin();
        }
    }

    T getValue() { return value; }
    void add(pair<int, int> d) {
        value.first += d.first;
        value.second += d.second;
    }
    int numChildren() { return children.size(); };
    int numRChildren() {
        unsigned int cnt = 0;
        for (unsigned int i = 0; i < children.size(); ++i) {
            cnt += children[i]->numRChildren();
        }
        return children.size() + cnt;
    };
    int importance() {
        unsigned int cnt = 0;
        for (unsigned int i = 0; i < children.size(); ++i) {
            cnt += children[i]->importance();
        }
        return value.third + cnt;
    };

    Node<T> *getChild(unsigned int i) {
        if (i < children.size()) { return children[i]; }
        else { return 0; }
    };

    void setChildren(vector< Node< T > * > new_children) {
        children = new_children;
    };
    void addChild(Node<T> *c) {
        children.push_back(c);
    }

    void addChild(T c) {
        Node<T> *n = new Node<T>(c);
        children.push_back(n);
    }

    void removeChild(unsigned int i) {
        if (i < children.size()) {
            children.erase(children.begin() + i);
        }
    }

    void removeRChild(unsigned int i) {
        if (i < children.size()) {
            typename vector<Node<T>* >::iterator n = children.begin();
            while (i--) { n++; }
            children.erase(n);
        }
    }

    Node<T> *operator[] (unsigned int i) {
        return getChild(i);
    }

};

#endif  /* NODE_HPP */

