
/**
 * @file
 * @author  Aapo Kyrola <akyrola@cs.cmu.edu>
 * @version 1.0
 *
 * @section LICENSE
 *
 * Copyright [2012] [Aapo Kyrola, Guy Blelloch, Carlos Guestrin / Carnegie Mellon University]
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 
 *
 * @section DESCRIPTION
 *
 * Vertex and Edge objects.
 */


#ifndef DEF_GRAPHCHI_OBJECTS
#define DEF_GRAPHCHI_OBJECTS

#include <vector>
#include <assert.h>
#include <omp.h>
#include <string.h>

#include "graphchi_types.hpp"
#include "util/qsort.hpp"

namespace graphchi {
    
/**
 * GNU COMPILER HACK TO PREVENT WARNINGS "Unused variable", if 
 * the particular app being compiled does not use a function.
 */
#ifdef __GNUC__
#define VARIABLE_IS_NOT_USED __attribute__ ((unused))
#else
#define VARIABLE_IS_NOT_USED
#endif
    
    
    
    template <typename EdgeDataType>
    class graphchi_edge {
        
    public:
        vid_t vertexid; // Source or Target vertex id. Clear from context.
        EdgeDataType * data_ptr;
        
        graphchi_edge() {}
        graphchi_edge(vid_t _vertexid, EdgeDataType * edata_ptr) : vertexid(_vertexid), data_ptr(edata_ptr) {
        }
        
#ifndef DYNAMICEDATA
        EdgeDataType get_data() {
            return * data_ptr;
        }
        
        void set_data(EdgeDataType x) {
            *data_ptr = x;
        }
#else 
        EdgeDataType * get_vector() {  // EdgeDataType is a chivector
            return data_ptr;
        }
#endif
        
        /**
          * Returns id of the endpoint of this edge. 
          */
        vid_t vertex_id() {
            return vertexid;
        }
        
 
    }  __attribute__((packed));
    
    template <typename ET>
    bool eptr_less(const graphchi_edge<ET> &a, const graphchi_edge<ET> &b) {
        return a.vertexid < b.vertexid;
    }
    
    
#ifdef SUPPORT_DELETIONS
    
    /*
     * Hacky support for edge deletions.
     * Edges are deleted by setting the value of the edge to a special
     * value that denotes it was deleted.
     * In the future, a better system could be designed.
     */
    
    // This is hacky...
    static inline bool VARIABLE_IS_NOT_USED is_deleted_edge_value(int val);
    static inline bool VARIABLE_IS_NOT_USED is_deleted_edge_value(bool val) {
        return val;
    }
    
    static inline bool VARIABLE_IS_NOT_USED is_deleted_edge_value(int val);
    static inline bool VARIABLE_IS_NOT_USED is_deleted_edge_value(int val) {
        assert(false); //cannot work
        return 0xffffffff == (unsigned int)val;
    }
    
    static inline bool VARIABLE_IS_NOT_USED is_deleted_edge_value(vid_t val);
    static inline bool VARIABLE_IS_NOT_USED is_deleted_edge_value(vid_t val) {
        return 0xffffffffffffffffu == val;
    }
    
    
    static inline bool VARIABLE_IS_NOT_USED is_deleted_edge_value(float val);
    static inline bool VARIABLE_IS_NOT_USED is_deleted_edge_value(float val) {
        return !(val < 0 || val > 0);
    }
    
    static void VARIABLE_IS_NOT_USED remove_edgev(graphchi_edge<bool> * e);
    static void VARIABLE_IS_NOT_USED remove_edgev(graphchi_edge<bool> * e) {
        e->set_data(true);
    }
    
    static void VARIABLE_IS_NOT_USED remove_edgev(graphchi_edge<vid_t> * e);
    static void VARIABLE_IS_NOT_USED remove_edgev(graphchi_edge<vid_t> * e) {
        e->set_data(0xffffffffffffffffu);
    }
    
    static void VARIABLE_IS_NOT_USED remove_edgev(graphchi_edge<int> * e);
    static void VARIABLE_IS_NOT_USED remove_edgev(graphchi_edge<int> * e) {
        e->set_data(0xffffffff);
    }
    
#endif  
    
    
    template <typename VertexDataType, typename EdgeDataType>
    class internal_graphchi_vertex {
        
    public:   // Todo, use friend
        // volatile int inc = 0;
        // volatile int outc = 0;
        
        vid_t vertexid = 0;

    protected:
        // graphchi_edge<EdgeDataType> * inedges_ptr = nullptr;
        // graphchi_edge<EdgeDataType> * outedges_ptr = nullptr;
        std::vector<graphchi_edge<EdgeDataType>> in_edges;
        std::vector<graphchi_edge<EdgeDataType>> out_edges;
        size_t nb_in_edges = 0; // sometimes we don't want to store the edges, just their number
        size_t nb_out_edges = 0;
        bool store_in_edges = true; // TODO use something more compact
        bool store_out_edges = true;

    public:
        bool modified = false;
        VertexDataType * dataptr = nullptr;


        /* Accessed directly by the engine */
        bool scheduled = false;
        bool parallel_safe = true;
        
#ifdef SUPPORT_DELETIONS
        int deleted_inc = 0;
        int deleted_outc = 0;
#endif
        
        
        internal_graphchi_vertex(): dataptr(nullptr){
#ifdef SUPPORT_DELETIONS
            deleted_outc = deleted_inc = 0;
#endif
        }
        
        internal_graphchi_vertex(vid_t _id) : vertexid(_id), store_in_edges(false), store_out_edges(false) {
        }
        
        internal_graphchi_vertex(vid_t _id,
                                 int indeg, 
                                 int outdeg) : vertexid(_id), nb_out_edges(outdeg) {
            in_edges.reserve(indeg);
            out_edges.reserve(outdeg);
        }

        internal_graphchi_vertex(vid_t _id,
                                 int indeg, 
                                 int outdeg, bool disable_outedge): 
                                 vertexid(_id), 
                                 nb_out_edges(outdeg), 
                                 store_out_edges(! disable_outedge) {
            in_edges.reserve(indeg);
            if (!disable_outedge) {
                out_edges.reserve(outdeg);
            }
        }
        
        virtual ~internal_graphchi_vertex() {}

        
        vid_t id() const {
            return vertexid;
        }
        
        size_t num_inedges() const { 
            if (store_in_edges) {
                return in_edges.size();
            }
            return nb_in_edges; 
            
        }
        size_t num_outedges() const {
            if (store_out_edges) {
                return out_edges.size();
            }
            return nb_out_edges; 
        }
        size_t num_edges() const { 
            return num_inedges() + num_outedges(); 
        }


        
                
        // Optimization: as only memshard (not streaming shard) creates inedgers,
        // we do not need atomic instructions here!
        inline void add_inedge(vid_t src, EdgeDataType * ptr, bool special_edge) {
            assert(src != vertexid);
#ifdef SUPPORT_DELETIONS
            // if (inedges_ptr != NULL && is_deleted_edge_value(*ptr)) {
            //     deleted_inc++;
            //     return;
            // }
#endif
            if (!store_in_edges) {
                nb_in_edges ++;
                return;
            }
            if (in_edges.capacity() <= in_edges.size()) {
                std::cout << "pb for vertex " << vertexid << " in capa = " << in_edges.capacity()  << " size " << in_edges.size() << std::endl;
            }
        //    assert(in_edges.capacity() < in_edges.size()); // we newer want to resize
            in_edges.push_back(graphchi_edge<EdgeDataType>(src, ptr));
            
          /*  if(inedges_ptr != NULL && inc > outedges_ptr - inedges_ptr) {
                logstream(LOG_FATAL) << "Tried to add more in-edges as the stored in-degree of this vertex (" << src << "). Perhaps a preprocessing step had failed?" << std::endl;
                assert(inc <= outedges_ptr - inedges_ptr);
            } */  // Deleted, since does not work when we have separate in-edge and out-edge arrays
        }
        
        inline void add_outedge(vid_t dst, EdgeDataType * ptr, bool special_edge) {
#ifdef SUPPORT_DELETIONS
            if (outedges_ptr != NULL && is_deleted_edge_value(*ptr)) {
                deleted_outc++;
                return;
            }
#endif
            if (!store_out_edges) {
                nb_out_edges ++;
                return;
            }
            if (out_edges.capacity() <= out_edges.size()) {
                std::cout << "pb for vertex " << vertexid << " out capa = " << out_edges.capacity()  << std::endl;
            }
        //    assert(out_edges.capacity() < out_edges.size()); // we newer want to resize
            out_edges.push_back(graphchi_edge<EdgeDataType>(dst, ptr));
            assert(dst != vertexid);
        }
        
        
    };
    
    template <typename VertexDataType, typename EdgeDataType >
    class graphchi_vertex : public internal_graphchi_vertex<VertexDataType, EdgeDataType> {
        
    public:
        
        graphchi_vertex() : internal_graphchi_vertex<VertexDataType, EdgeDataType>() { }
        
        graphchi_vertex(vid_t _id,
                        int indeg, 
                        int outdeg,
                        bool disable_out_edges) : 
            internal_graphchi_vertex<VertexDataType, EdgeDataType>(_id, indeg, outdeg, disable_out_edges) {}

        graphchi_vertex(vid_t _id) : 
            internal_graphchi_vertex<VertexDataType, EdgeDataType>(_id) {}

        virtual ~graphchi_vertex() {}
        
        /** 
          * Returns ith edge of a vertex, ignoring 
          * edge direction.
          */
        graphchi_edge<EdgeDataType> * edge(size_t i) {
            if (i < this->inc) return inedge(i);
            else return outedge(i - this->inc);
        }

        
        graphchi_edge<EdgeDataType> * inedge(size_t i) {
            assert(i >= 0 && i < this->in_edges.size());
            return &this->in_edges[i];
        }
        
        graphchi_edge<EdgeDataType> * outedge(size_t i) {
            assert(i >= 0 && i < this->out_edges.size());
            return &this->out_edges[i];
        }        
        
        graphchi_edge<EdgeDataType> * random_outedge() {
            if (this->out_edges.empty()) return NULL;
            return outedge((std::abs(random()) % this->out_edges.size()));
        }
            
        /** 
          * Get the value of vertex
          */
#ifndef DYNAMICVERTEXDATA
        VertexDataType get_data() {
            return *(this->dataptr);
        }
#else
        // VertexDataType must be a chivector
        VertexDataType * get_vector() {
            this->modified = true;  // Assume vector always modified... Temporaryh solution.
            return this->dataptr;
        }
#endif
        
        /**
          * Modify the vertex value. The new value will be
          * stored on disk.
          */
        virtual void set_data(VertexDataType d) {
            *(this->dataptr) = d;
            this->modified = true;
        }
        
        // TODO: rethink
        static bool computational_edges() {
            return false;
        }
        static bool read_outedges() {
            return true;
        }
        
        
        /**
         * Sorts all the edges. Note: this will destroy information
         * about the in/out direction of an edge. Do use only if you
         * ignore the edge direction.
         */
        // void VARIABLE_IS_NOT_USED sort_edges_indirect() {
        //     // Check for deleted edges first...
        //     if (this->inc != (this->outedges_ptr - this->inedges_ptr)) {
        //         // Moving
        //         memmove(&this->inedges_ptr[this->inc], this->outedges_ptr, this->outc * sizeof(graphchi_edge<EdgeDataType>));
        //         this->outedges_ptr = &this->inedges_ptr[this->inc];
        //     }
        //     quickSort(this->inedges_ptr, (int) (this->inc + this->outc), eptr_less<EdgeDataType>);
            
        // }
        
        
#ifdef SUPPORT_DELETIONS
        void VARIABLE_IS_NOT_USED remove_edge(int i) {
            do not compil
            remove_edgev(edge(i));
        }
        
        void VARIABLE_IS_NOT_USED remove_inedge(int i) {
            remove_edgev(inedge(i));
        }
        
        void VARIABLE_IS_NOT_USED remove_outedge(int i) {
            remove_edgev(outedge(i));
        }
        
        void VARIABLE_IS_NOT_USED remove_alledges() {
            for(int j=this->num_edges()-1; j>=0; j--) remove_edge(j);
        }
        
        
        void VARIABLE_IS_NOT_USED remove_alloutedges() {
            for(int j=this->num_outedges()-1; j>=0; j--) remove_outedge(j);
        }
        
        void VARIABLE_IS_NOT_USED remove_allinedges() {
            for(int j=this->num_inedges()-1; j>=0; j--) remove_inedge(j);
        }
        
        
#endif
        
        
    };
    
    /**
      * Experimental code
      */
    
    // If highest order bit is set, the edge is "special". This is used
    // to indicate - in the neighborhood model - that neighbor's value is
    // cached in memory. 
#define HIGHMASK (1 + (2147483647 >> 1))
#define CLEARMASK (2147483647 >> 1)
    inline vid_t translate_edge(vid_t rawid, bool &is_special) {
        is_special = (rawid & HIGHMASK) != 0;
        return rawid & CLEARMASK;
    }
    inline vid_t make_special(vid_t rawid) {
        return rawid | HIGHMASK;
    }
    inline bool is_special(vid_t rawid) {
        return (rawid & HIGHMASK) != 0;
    }
    
    

} // Namespace

#endif
