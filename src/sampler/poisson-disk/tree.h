// $Id: tree.h 790 2011-01-06 19:10:37Z mag $
#ifndef TREE_H_
#define TREE_H_

#include <stack>
#include <cmath>
#include <vector>
#include <cassert>
#include <boost/pool/pool.hpp>

#include "cube.h"
#include "tuple.h"
#include "sample.h"
#include "bignum.h"
#include "interval.h"

using namespace boost;

template <unsigned char N>
class Tree
{

    struct node
    {

        static unsigned char lmax;

        virtual ~node() {}

        virtual BigNum<N> volume(unsigned char l) const = 0;
        virtual bool valid(unsigned char l) const = 0;
        virtual bool leaf() const = 0;

    };

    struct inner_node : public node
    {

        static pool<> mempool;
        static uint32_t allocs;

        node* _c[2];
        BigNum<N> _v;

        inner_node();
        inner_node(const BigNum<N>& v,
                   node* c[2]);
        virtual ~inner_node() {}

        void prune();

        BigNum<N> update(const Cube<N>& c,
                         const Sample<N>& s,
                         uint32_t i,
                         unsigned char l);
        BigNum<N> generate(const Cube<N>& c,
                           Sample<N>& sample,
                           unsigned char l);

        void push(stack<const inner_node*>& s) const;

        virtual BigNum<N> volume(unsigned char l) const { return _v; }
        virtual bool valid(unsigned char l) const;
        virtual bool leaf() const { return false; }

        void* operator new (size_t);
        void operator delete (void* b);

    };

    struct leaf_node : public node
    {

        static bool sph;

        static BigNum<N>* vol;
        static Cube<N>* bound;

        static const unsigned char DMAX = 24;
        static const unsigned char LMAX = DMAX*N;

        static pool<> mempool;
        static uint32_t allocs;
        static uint32_t terminals;

        struct sample
        {

            static pool<> mempool;
            static uint32_t allocs;

            uint32_t _i;
            const sample* _n;

            sample(uint32_t i,
                   const sample* n);

            void* operator new (size_t);
            void operator delete (void* b);

        };

        const sample* _s;

        leaf_node();
        leaf_node(const sample* s);
        virtual ~leaf_node();

        void append(uint32_t i);

        BigNum<N> subdivide(const Cube<N>& c,
                            unsigned char l,
                            node* n[2]);

        BigNum<N> update(const Cube<N>& c,
                         const Sample<N>& s,
                         uint32_t i,
                         unsigned char l,
                         node* n[2]);
        BigNum<N> generate(const Cube<N>& c,
                           Sample<N>& sample,
                           unsigned char l,
                           node* n[2]);

        bool valid(const Sample<N>& s) const;

        virtual BigNum<N> volume(unsigned char l) const { return vol[l]; }
        virtual bool valid(unsigned char) const { return true; }
        virtual bool leaf() const { return true; }

        void* operator new (size_t);
        void operator delete (void* b);

    };

    node* _r;

    Cube<N> _c;

    struct NodeContext
    {
       node* m_node;
       Cube<N> m_cube;
       unsigned char m_level;
       NodeContext() :
       m_node(0), m_level(0UL) {}
       NodeContext(node* node,
                   const Cube<N>& cube,
                   unsigned char level) :
       m_node(node), m_cube(cube), m_level(level) {}
       friend bool operator == (const NodeContext& n1,
                                const NodeContext& n2)
       { return n1.m_level == n2.m_level &&
                n1.m_cube == n2.m_cube &&
                n1.m_node == n2.m_node; }
    };

    typedef stack<NodeContext> NodeContextStack;

public:

    class Iterator: private NodeContextStack
    {

        Iterator(node* n,
                 const Cube<N>& c,
                 unsigned char l);
        Iterator(const NodeContextStack& s);

        void next();

        friend class Tree;

    public:

        Iterator();

        Cube<N> operator * () const;

        Iterator& operator ++ (int);

        bool operator != (const Iterator& i) const;

    };

    Tree();
    ~Tree();

    Sample<N> generate();

    bool update(const Sample<N>& s,
                uint32_t i);

    bool empty() const;

    bool valid() const;

    double volume() const;

    Iterator end() const;
    Iterator begin() const;

    uint32_t get_terminals() const;

    void spherical_boundary();
    void box_boundary(const Tuple<N,float>& bound);

    static void memory_stats();

};

template <unsigned char N>
unsigned char
Tree<N>::node::lmax = 0U;

template <unsigned char N>
boost::pool<>
Tree<N>::inner_node::mempool(sizeof(inner_node));

template <unsigned char N>
uint32_t
Tree<N>::inner_node::allocs = 0UL;

template <unsigned char N>
bool
Tree<N>::leaf_node::sph = false;

template <unsigned char N>
uint32_t
Tree<N>::leaf_node::terminals = 0UL;

template <unsigned char N>
Cube<N>*
Tree<N>::leaf_node::bound = 0;

template <unsigned char N>
boost::pool<>
Tree<N>::leaf_node::mempool(sizeof(leaf_node));

template <unsigned char N>
uint32_t
Tree<N>::leaf_node::allocs = 0UL;

template <unsigned char N>
boost::pool<>
Tree<N>::leaf_node::sample::mempool(sizeof(sample));

template <unsigned char N>
uint32_t
Tree<N>::leaf_node::sample::allocs = 0UL;

template <unsigned char N>
BigNum<N>*
Tree<N>::leaf_node::vol = 0;

template <unsigned char N>
inline bool
Tree<N>::empty() const
{ return !_r; }

template <unsigned char N>
inline typename Tree<N>::Iterator
Tree<N>::begin() const
{ return Iterator(_r, _c, 0); }

template <unsigned char N>
inline typename Tree<N>::Iterator
Tree<N>::end() const
{ return Iterator(); }

template <unsigned char N>
Tree<N>::inner_node::inner_node()
{ _c[0] = _c[1] = 0; }

template <unsigned char N>
Tree<N>::leaf_node::leaf_node() : _s(0) {}

template <unsigned char N>
inline
Tree<N>::inner_node::inner_node(const BigNum<N>& v,
                                node* c[2]) :
    _v(v) { _c[0] = c[0]; _c[1] = c[1]; }

template <unsigned char N>
inline
Tree<N>::leaf_node::leaf_node(const sample* s) :
    _s(s) {}

template <unsigned char N>
inline void
Tree<N>::leaf_node::append(uint32_t i)
{ _s = new sample(i,_s); }

template <unsigned char N>
void
Tree<N>::inner_node::push(stack<const inner_node*>& s) const
{

   for (unsigned char i = 0; i < 2; ++i)
    {
        if (_c[i])
        {
           if (!_c[i]->leaf())
              s.push(reinterpret_cast<const inner_node*>(_c[i]));
           else
              delete _c[i];
        }
    }

}

template <unsigned char N>
inline
Tree<N>::Iterator::Iterator() {}

template <unsigned char N>
inline
Tree<N>::Iterator::Iterator(node* n,
                            const Cube<N>& c,
                            unsigned char l)
{ if (n) { NodeContextStack::push(NodeContext(n, c, l)); next(); } }

template <unsigned char N>
inline
Tree<N>::Iterator::Iterator(const NodeContextStack& s) :
    NodeContextStack(s) {}

template <unsigned char N>
inline Cube<N>
Tree<N>::Iterator::operator * () const
{ return !NodeContextStack::empty() ? NodeContextStack::top().m_cube : Cube<N>(); }

template <unsigned char N>
inline typename Tree<N>::Iterator&
Tree<N>::Iterator::operator ++ (int)
{ if (!NodeContextStack::empty()) { NodeContextStack::pop(); next(); }  return *this; }

template <unsigned char N>
inline bool
Tree<N>::Iterator::operator != (const Iterator& i) const
{ const NodeContextStack& it = reinterpret_cast<const NodeContextStack&>(i);
  return it != *this; }

template <unsigned char N>
void
Tree<N>::Iterator::next()
{

    while (!NodeContextStack::empty())
    {

        NodeContext nc(NodeContextStack::top());

        if (nc.m_node->leaf())
            return;

        NodeContextStack::pop();

        inner_node* n = reinterpret_cast<inner_node*>(nc.m_node);

        for (unsigned char i = 0; i < 2; ++i)
        {
            if (n->_c[i])
               push(NodeContext(n->_c[i], nc.m_cube.child(nc.m_level, i), nc.m_level + 1));
        }

    }

}

template <unsigned char N>
inline
Tree<N>::leaf_node::sample::sample(uint32_t i,
                                   const sample* n) :
    _i(i), _n(n) {}

template <unsigned char N>
Tree<N>::Tree() :
    _r(new leaf_node()),
    _c(Tuple<N,Interval>(Interval()))
{

   leaf_node::vol = new BigNum<N>[leaf_node::LMAX + 1];

   for (unsigned char l = leaf_node::LMAX + 1; l > 0; l--)
   {

      unsigned char i = l - 1;

      if (i % N)
      {
         assert(i+1 <= leaf_node::LMAX);
         BigNum<N> n(leaf_node::vol[i+1]);
         n += n; leaf_node::vol[i] = n;
      }
      else
         leaf_node::vol[i] = BigNum<N>(leaf_node::DMAX - i/N);

   }

   node::lmax = static_cast<unsigned char>(-logf(Sample<N>::getRadius()/sqrtf(N))/logf(2.0f)*N);

}

template <unsigned char N>
Tree<N>::~Tree()
{
   delete [] leaf_node::vol;
   delete leaf_node::bound;
}

template <unsigned char N>
void
Tree<N>::spherical_boundary()
{
   leaf_node::sph = true;
}

template <unsigned char N>
void
Tree<N>::box_boundary(const Tuple<N,float>& bound)
{

   leaf_node::bound = new Cube<N>();

   for (unsigned char i = 0; i < N; ++i)
   {
      double size = bound.get(i);
      leaf_node::bound->get(i) = Interval(0.5*(1.0 - size), 0.5*(1.0 + size));
   }

}

template <unsigned char N>
double
Tree<N>::volume() const
{

    if (_r == 0)
        return 1.0;

    if (_r->leaf())
        return 0.0;

    inner_node* n = reinterpret_cast<inner_node*>(_r);

    return 1.0 - static_cast<double>(n->_v);

}

template <unsigned char N>
Sample<N>
Tree<N>::generate()
{

    Sample<N> sample;

    if (_r)
    {

        if (_r->leaf())
        {

            leaf_node* l = reinterpret_cast<leaf_node*>(_r);

            Sample<N> s(_c.randomise());

            if (!l->valid(s))
            {

                node* n[2];
                BigNum<N> lv(l->volume(0));
                BigNum<N> dv(l->subdivide(_c, 0, n));
            
                assert(!lv || dv <= lv);

                if (lv && dv == lv)
                {
                    delete l;
                    _r = 0;
                }
                else
                {
                    lv -= dv;
                    _r = new inner_node(lv, n);
                    delete l;
                }

            }
            else
                sample = s;
 
        }
        else
        {

            inner_node* n = reinterpret_cast<inner_node*>(_r);
            BigNum<N> dv(n->generate(_c, sample, 0));
            assert(!n->_v || dv <= n->_v);

            if (n->_v && dv == n->_v)
            {
                delete n;
                _r = 0;
            }
            else
                n->_v -= dv;

        }

    }

    return sample;

}

template <unsigned char N>
BigNum<N>
Tree<N>::inner_node::generate(const Cube<N>& c, Sample<N>& sample, unsigned char l)
{

    unsigned char l1 = l + 1;
    BigNum<N> pv(this->_v.random());
    unsigned char i_next = (_c[0] && _c[0]->volume(l1) > pv) ? 0 : 1;
    assert(_c[i_next]);

    if (_c[i_next]->leaf())
    {

        node* n[2];
        leaf_node* l_next = reinterpret_cast<leaf_node*>(_c[i_next]);
        BigNum<N> lv(l_next->volume(l1));
        BigNum<N> dv(l_next->generate(c.child(l, i_next), sample, l1, n));
        assert(dv <= lv);

        if (dv < lv)
        {
           lv -= dv;
           if (n[0] || n[1])
           {
              _c[i_next] = new inner_node(lv, n);
              delete l_next;
           }
        }
        else
        {
           _c[i_next] = 0;        
           delete l_next;
        }

        return dv;

    }

    inner_node* n_next = reinterpret_cast<inner_node*>(_c[i_next]);
    BigNum<N> dv(n_next->generate(c.child(l, i_next), sample, l1));
    assert(dv <= n_next->_v);

    if (n_next->_v == dv)
    {
       _c[i_next] = 0;
       delete n_next;
    }
    else
       n_next->_v -= dv;

    return dv;

}

template <unsigned char N>
BigNum<N>
Tree<N>::leaf_node::generate(const Cube<N>& c, Sample<N>& sample, unsigned char l, node* n[2])
{

    Sample<N> s(c.randomise());

    if (valid(s))
    {
       sample = s;
       n[0] = n[1] = 0;
       return BigNum<N>();
    }

    if (l < LMAX)
       return subdivide(c, l, n);

    ++terminals;

    return this->volume(l);

}

template <unsigned char N>
void
Tree<N>::inner_node::prune()
{

    stack<const inner_node*> nstack;

    push(nstack);
    while (!nstack.empty()) {
        const inner_node* n = nstack.top();
        nstack.pop();
        n->push(nstack);
        delete n;
    }

    _c[0] = _c[1] = 0;

}

template <unsigned char N>
bool
Tree<N>::leaf_node::valid(const Sample<N>& s) const
{

    if (bound)
    {

       if (!bound->contains(s))
          return false;

    }

    if (sph)
    {
       
       double d = 0.0;

       for (unsigned char i = 0; i < N; ++i)
       {

           double di = s.get(i) - 0.5;

           d += di*di;
           
           if (d > 0.25)
              return false;

       }

    }

    if (_s)
    {
        const sample* ns = _s;
        while (ns)
        {
            if (Sample<N>::getSample(ns->_i).contains(s))
                return false;
            ns = ns->_n;
        }
    }

    return true;

}

template <unsigned char N>
bool
Tree<N>::update(const Sample<N>& s, uint32_t i)
{

    typename Sample<N>::Status status = s.intersects(_c);

    if (status == Sample<N>::Out)
        return false;

    if (_r)
    {

        if (_r->leaf())
        {

            leaf_node* l = reinterpret_cast<leaf_node*>(_r);

            if (status == Sample<N>::In)
            {
                delete l;
                _r = 0;
            }
            else
            {

                l->append(i);

                node* n[2];
                BigNum<N> lv(l->volume(0));
                BigNum<N> dv(l->subdivide(_c, 0, n));

                assert(!lv || dv <= lv);

                if (lv && dv == lv)
                {
                    assert(n[0] == 0 && n[1] == 0);
                    delete l;
                    _r = 0;
                }
                else
                {
                    lv -= dv;
                    _r = new inner_node(lv, n);
                    delete l;
                }

            }

        }
        else
        {

            inner_node* n = reinterpret_cast<inner_node*>(_r);
            BigNum<N> dv(n->update(_c, s, i, 0));
            assert(!n->_v || dv <= n->_v);

            if (n->_v && dv == n->_v)
            {
                delete n;
                _r = 0;
            }
            else
                n->_v -= dv;

        }

    }

    return true;

}

template <unsigned char N>
BigNum<N>
Tree<N>::inner_node::update(const Cube<N>& c, const Sample<N>& s, uint32_t i, unsigned char l)
{

    if (l > node::lmax && c.contains(s))
    {
       prune();
       return this->_v;
    }

    typename Sample<N>::Status status = s.intersects(c);

    if (status == Sample<N>::Out)
        return BigNum<N>();

    if (status == Sample<N>::In)
    {
        prune();
        return this->_v;
    }

    BigNum<N> dv;
    Cube<N> cube[2];

    c.subdivide(l, cube);

    l += 1;

    for (unsigned char j = 0; j < 2; ++j)
    {

       if (_c[j])
       {

          if (_c[j]->leaf())
          {

             node* n[2];
             leaf_node* l_curr = reinterpret_cast<leaf_node*>(_c[j]);
             BigNum<N> lv(l_curr->volume(l));
             BigNum<N> dvl(l_curr->update(cube[j], s, i, l, n));
             assert(dvl <= lv);

             if (dvl < lv)
             {
                lv -= dvl;
                if (n[0] || n[1])
                {
                   _c[j] = new inner_node(lv, n);
                   delete l_curr;
                }
             }
             else
             {
                _c[j] = 0;
                delete l_curr;
             }

             dv += dvl;

          }
          else
          {

             inner_node* n_curr = reinterpret_cast<inner_node*>(_c[j]);
             BigNum<N> dvc(n_curr->update(cube[j], s, i, l));
             assert(dvc <= n_curr->_v);

             if (n_curr->_v == dvc)
             {
                _c[j] = 0;
                delete n_curr;
             }
             else
                n_curr->_v -= dvc;

             dv += dvc;

          }

       }

    }

    return dv;

}

template <unsigned char N>
BigNum<N>
Tree<N>::leaf_node::update(const Cube<N>& c, const Sample<N>& s, uint32_t i, unsigned char l, node* n[2])
{

    if (l > node::lmax && c.contains(s))
        return this->volume(l);

    typename Sample<N>::Status status = s.intersects(c);

    if (status == Sample<N>::Out)
    {
        n[0] = n[1] = 0;
        return BigNum<N>();
    }

    if (status == Sample<N>::In)
       return this->volume(l);

    append(i);

    return subdivide(c, l, n);

}

template <unsigned char N>
inline bool
Tree<N>::valid() const
{

    if (_r == 0)
        return inner_node::allocs == 0UL &&
               leaf_node::allocs == 0UL &&
               leaf_node::sample::allocs == 0UL;

    if (_r->leaf())
        return _r->volume(0) == BigNum<N>();

    return _r->valid(0);

}

template <unsigned char N>
bool
Tree<N>::inner_node::valid(unsigned char l) const
{

    l += 1;

    BigNum<N> v;

    for (unsigned char i = 0; i < 2; ++i)
    {
       if (_c[i])
         v += _c[i]->volume(l);
    }

    if (this->_v != v)
        return false;

    for (unsigned char i = 0; i < 2; ++i)
    {
       if (_c[i] && !_c[i]->valid(l))
          return false;
    }

    return true;

}

template <unsigned char N>
Tree<N>::leaf_node::~leaf_node()
{

    while (_s)
    {
        const sample* s = _s;
        _s = s->_n;
        delete s;
    }

}

template <unsigned char N>
BigNum<N>
Tree<N>::leaf_node::subdivide(const Cube<N>& c, unsigned char l, node* n[2])
{

    Cube<N> cube[2];

    c.subdivide(l, cube);

    BigNum<N> dv;
    BigNum<N> dvc(this->volume(l+1));

    for (unsigned char i = 0; i < 2; i++)
    {
        
        if (bound)
        {

            if (!bound->overlaps(cube[i]))
            {
               dv += dvc;
               n[i] = 0;
               continue;
            }

        }

        if (sph)
        {

           double d = 0.0;

           for (unsigned char j = 0; j < N; ++j)
           {

              double min = cube[i].get(j).min() - 0.5;

              if (min > 0.0) {

                 d += min*min;

                 if (d > 0.25)
                    break;

                 continue;

              }

              min = 0.5 - cube[i].get(j).max();

              if (min > 0.0)
              {

                  d += min*min;

                  if (d > 0.25)
                     break;

                  continue;

              }

           }

           if (d > 0.25)
           {
               dv += dvc;
               n[i] = 0;
               continue;
           }

        }

        const sample* s = _s;
        const sample* sc = 0;

        while (s) {

            typename Sample<N>::Status status = 
                Sample<N>::getSample(s->_i).intersects(cube[i]);

            if (status == Sample<N>::In)
                break;

            if (status == Sample<N>::Over)
                sc = new sample(s->_i, sc);

            s = s->_n;

        }

        if (s) {

            while (sc) {
                const sample* s = sc;
                sc = s->_n;
                delete s;
            }

            dv += dvc;
            n[i] = 0;
            continue;

        }

        n[i] = new leaf_node(sc);

    }

    while (_s) {
        const sample* s = _s;
        _s = s->_n;
        delete s;
    }

    return dv;

}

template <unsigned char N>
void*
Tree<N>::inner_node::operator new (size_t)
{

    void* p = 0;

    while (p == 0)
    {

        p = mempool.malloc();

        if (p == 0)
        {

            pool<>::size_type next_size = mempool.get_next_size() >> 1;

            if (next_size < sizeof(node))
            {
                cerr << "\nTree<N>::inner_node: Unable to allocate more memory. Exiting...\n";
                exit(EXIT_FAILURE);
            }
            else
                mempool.set_next_size(next_size);

        }

    }

    ++allocs;

    return p;

}

template <unsigned char N>
void*
Tree<N>::leaf_node::operator new (size_t)
{

    void* p = 0;

    while (p == 0)
    {

        p = mempool.malloc();

        if (p == 0)
        {

            pool<>::size_type next_size = mempool.get_next_size() >> 1;

            if (next_size < sizeof(leaf_node))
            {
                cerr << "\nTree<N>::leaf_node: Unable to allocate more memory. Exiting...\n";
                exit(EXIT_FAILURE);
            }
            else
                mempool.set_next_size(next_size);

        }

    }

    ++allocs;

    return p;

}

template <unsigned char N>
inline void
Tree<N>::inner_node::operator delete (void* b)
{ mempool.free(b); --allocs; }

template <unsigned char N>
inline void
Tree<N>::leaf_node::operator delete (void* b)
{ mempool.free(b); --allocs; }

template <unsigned char N>
inline void*
Tree<N>::leaf_node::sample::operator new (size_t)
{

    void* p = 0;

    while (p == 0)
    {

        p = mempool.malloc();

        if (p == 0)
        {

            pool<>::size_type next_size = mempool.get_next_size() >> 1;

            if (next_size < sizeof(typename leaf_node::sample))
            {
                cerr << "\nTree<N>::leaf_node::sample: Unable to allocate more memory. Exiting...\n";
                exit(EXIT_FAILURE);
            }
            else
                mempool.set_next_size(next_size);

        }

    }

    ++allocs;

    return p;

}

template <unsigned char N>
inline void
Tree<N>::leaf_node::sample::operator delete (void* b)
{ mempool.free(b); --allocs; }

template <unsigned char N>
inline uint32_t
Tree<N>::get_terminals() const
{ return leaf_node::terminals; }

template <unsigned char N>
inline void
Tree<N>::memory_stats()
{ cout << inner_node::allocs << ' ' << leaf_node::allocs << ' ' << leaf_node::sample::allocs << endl; }

#endif /*TREE_H_*/
