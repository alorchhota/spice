#include <Rcpp.h>
using namespace Rcpp;

NumericVector get_lower_triangle_vector(const NumericMatrix x, const bool diagonal, const bool bycol);

size_t find_component(size_t node,
                      std::vector<size_t> &node_2_component,
                      std::map<size_t, std::vector<size_t> > &component_2_nodes,
                      size_t &last_component){
  size_t comp = node_2_component.at(node);
  // create if does not exist
  if(comp == 0){
    last_component++;
    node_2_component.at(node) = last_component;
    std::vector<size_t> nodes(1,node);
    component_2_nodes[last_component] = nodes;
    comp = last_component;
  }
  return comp;
}

void union_component(size_t node1,
                     size_t node2,
                     std::vector<size_t> &node_2_component,
                     std::map<size_t, std::vector<size_t> > &component_2_nodes){
  size_t c1 = node_2_component.at(node1);
  size_t c2 = node_2_component.at(node2);
  std::vector<size_t>& c1_nodes = component_2_nodes[c1];
  std::vector<size_t>& c2_nodes = component_2_nodes[c2];
  // set c2 nodes to c1 component
  for(size_t i=0; i<c2_nodes.size(); i++){
    size_t c2_node = c2_nodes.at(i);
    node_2_component.at(c2_node) = c1;
  }

  // update c1 component nodes
  c1_nodes.insert(c1_nodes.end(), c2_nodes.begin(), c2_nodes.end() );

  // delete c2 component nodes
  c2_nodes.resize(0);
}

template<typename T>
std::pair<T, T> lowerTriangleVectorIdx2matrixIdx(T idx){
  T i = ceil(sqrt(2.0 * idx + 2.25) - 0.5);
  T j = idx - (i - 1.0) * i / 2;
  return std::make_pair(i, j);
}


// [[Rcpp::export]]
DataFrame kruskal_mst_c(const NumericMatrix x, const bool maximum=false, const bool verbose = false){
  size_t n_gene = x.nrow(), nc = x.ncol();
  if(n_gene != nc)
    stop("x must be a square matrix.");

  std::vector<size_t> node_2_component(n_gene, 0);
  std::map<size_t, std::vector<size_t> > component_2_nodes;
  size_t last_component = 0;
  std::vector<double> MST_from(n_gene-1);
  std::vector<double> MST_to(n_gene-1);
  std::vector<double> MST_weight(n_gene-1);

  NumericVector weights = get_lower_triangle_vector(x, false, false);
  Function rorder("order");
  NumericVector weight_orders = rorder(weights, _["decreasing"] = maximum, _["na.last"] = NA_LOGICAL);

  size_t n_edges = 0;
  size_t si = 0;
  while(si < weight_orders.size()){
    size_t weight_idx = weight_orders.at(si) -1 ; // -1 for 0-based index in c++
    std::pair<size_t, size_t> mat_idx = lowerTriangleVectorIdx2matrixIdx(weight_idx);
    size_t node1 = mat_idx.first;
    size_t node2 = mat_idx.second;
    size_t c1 = find_component(node1, node_2_component, component_2_nodes, last_component);
    size_t c2 = find_component(node2, node_2_component, component_2_nodes, last_component);
    if (c1 != c2){
      MST_from.at(n_edges) = node1+1;  // +1 for 1-based index in R
      MST_to.at(n_edges) = node2+1;    // +1 for 1-based index in R
      MST_weight.at(n_edges) = x(node1, node2);
      union_component(node1, node2, node_2_component, component_2_nodes);
      n_edges = n_edges + 1;
      if(verbose && n_edges == n_gene-1){
        Rcout << "MST completed at rank " << si+1 << " (max " << weight_orders.size() << ")." << std::endl;
        break;
      }
    }
    si++;
  }

  if(n_edges < n_gene-1){
    MST_from.resize(n_edges);
    MST_to.resize(n_edges);
    MST_weight.resize(n_edges);
  }

  DataFrame MST_df = DataFrame::create(_["from"] = MST_from,
                                       _["to"] = MST_to,
                                       _["weight"] = MST_weight,
                                       _["stringsAsFactors"] = false);
  return MST_df;
}
