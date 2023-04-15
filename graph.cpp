/**
 * @brief A graph is made up of vertices and edges.
 * Vertex labels are unique. A vertex can be connected to other vertices via
 * weighted, directed edges or nondirected egdes. A vertex cannot
 * connect to itself or have multiple edges to the same vertex. Can traverse
 * graph using bfs or dfs, find shortest path between vertices with Dijkstra
 * algorithm.
 */

#include "graph.h"
#include <algorithm>
#include <climits>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>
using namespace std;

/**
 * constructor, empty map "graph" & bool "directionalEdges"
 * map "graph" has key-value pair <string label, vector<string> adjVertices>
 * "directionalEdges" defaults as true
 * @param directionalEdges bool to determine whether edges are directional
 */
Graph::Graph(bool directionalEdges) {
  this->directionalEdges = directionalEdges;
}

/**
 * destructor
 */
Graph::~Graph() {}

/**
 * each unique vertex has its own entry in the graph
 * @return int total number of vertices
 */
int Graph::verticesSize() const {
  int size = theGraph.size();
  return size;
}

/**
 * only single edges between vertices so just traverse entries
 * and count number of outgoing edges (size of adjVertices vector>
 * @return int total number of edges in a directed graph
 */
int Graph::directionalEdgesCounter() const {
  int total = 0;
  for (auto entries : theGraph) {
    total += entries.second.size();
  }
  return total;
}

/**
 * edges point both directions in an undirected graph, so count
 * number of outgoing edges for each vertice (size of adjVertice vector)
 * then divide by 2
 * @return int total number of edges in undirected graph
 */
int Graph::nondirectionalEdgesCounter() const {
  int total = 0;
  for (auto entries : theGraph) {
    total += entries.second.size();
  }
  return (total / 2);
}

/**
 * if the graph is empty there is no edges, else call the
 * corresponding helper function for the type of graph
 * @return int total number of edges in (non)directional graph
 */
int Graph::edgesSize() const {
  if (theGraph.size() == 0) {
    return 0;
  }

  if (directionalEdges) {
    return directionalEdgesCounter();
  } else {
    return nondirectionalEdgesCounter();
  }
}

/**
 * the adjVector (value in key-value pair) contains all outgoing edges
 * of a vertex, so return size of adjVector
 * if the vertex is not found, -1 is returned
 * @param label: the vertex being searched for
 * @return int total number of outgoing edges from a vertex
 */
int Graph::vertexDegree(const string &label) const {
  if (contains(label)) {
    return theGraph.at(label).size();
  } else {
    return -1;
  }
}

/**
 * add the key-value pair to the map, if the key (label) is already
 * in the map return false
 * @param label: the vertex being searched for
 * @return bool whether or not the vertex was added
 */
bool Graph::add(const string &label) {
  if (!contains(label)) {
    vector<vector<string>> friends;
    theGraph.insert({label, friends});
    return true;
  } else {
    return false;
  }
}

/**
 * uses the "find" function of std::map to search for label
 * if the matching label is not found return false
 * @param label: the vertex being searched for
 * @return bool whether or not the vertex is in the graph
 */
bool Graph::contains(const string &label) const {
  if (theGraph.find(label) == theGraph.end()) {
    return false;
  } else {
    return true;
  }
}

/**
 * traverse the adjVector of the desired vertex (the value of the key-value
 * pair) and add label(weight) string to a vector. sort the vector then
 * concatenate all elements to the ret string
 * to a string (A-3->B, A-5->C should return B(3),C(5))
 * @param label: the vertex being searched for
 * @return string of adjacent vertices and their weights
 */
string Graph::getEdgesAsString(const string &label) const {
  string answer = "";
  vector<string> answerStorage;
  for (auto pair : theGraph.at(label)) {
    string weight = pair[0];
    string destination = pair[1];
    string currSection = ((destination) + "(" + (weight) + "),");
    answerStorage.push_back(currSection);
  }
  sort(answerStorage.begin(), answerStorage.end());
  for (auto elm : answerStorage) {
    answer += elm;
  }
  return answer.substr(0, answer.length() - 1);
}

/**
 * add the desired edge only to the "from" vertex adjVector,
 * since directed graphs can only have one edge between two vertices
 * @param from: the starting vertex
 * @param to: the destination vertex
 * @param weight: the weight of the edge between the two vertices
 * @return bool whether or not the edge was successfully added
 */
bool Graph::directionalConnection(const string &from, const string &to,
                                  int weight) {
  // check if there is an existing edge between the two
  // going the opposite direction
  for (auto pair : theGraph.at(to)) {
    if (pair[1] == from) {
      return false;
    }
  }
  // check if there is an existing edge between the two
  // going the same direction
  for (auto pair : theGraph.at(from)) {
    if (pair[1] == to) {
      return false;
    }
  }
  // there's no matching edges, so we can add a new one
  // going the specified direction
  vector<string> newConnection;
  newConnection.push_back(to_string(weight));
  newConnection.push_back(to);
  theGraph[from].push_back(newConnection);
  // sort the connections
  sort(theGraph[from].begin(), theGraph[from].end());
  return true;
}

/**
 * add the desired edge to both the "from" and "to" vertices
 * since an undirected graph's edges must point both ways
 * @param from: the starting vertex
 * @param to: the destination vertex
 * @param weight: the weight of the edge between the two vertices
 * @return bool whether or not the edge was successfully added
 */
bool Graph::nondirectionalConnection(const string &from, const string &to,
                                     int weight) {

  // check if there is an existing edge
  for (auto pair : theGraph.at(from)) {
    if (pair[1] == to) {
      return false;
    }
  }
  // there's no matching edges, so we can add a new one
  // going both directions
  // from -> to
  vector<string> newConnectionFrom;
  newConnectionFrom.push_back(to_string(weight));
  newConnectionFrom.push_back(to);
  theGraph[from].push_back(newConnectionFrom);
  // sort the connections
  sort(theGraph[from].begin(), theGraph[from].end());
  // to -> from
  vector<string> newConnectionTo;
  newConnectionTo.push_back(to_string(weight));
  newConnectionTo.push_back(from);
  theGraph[to].push_back(newConnectionTo);
  // sort the connections
  sort(theGraph[to].begin(), theGraph[to].end());
  return true;
}

/**
 * cannot create a self edge, add the vertex if it is missing
 * call the corresponding helper function for (non)directional graph
 * @param from: the starting vertex
 * @param to: the destination vertex
 * @param weight: the weight of the edge between the two vertices
 * @return bool whether or not the edge was successfully added
 */
bool Graph::connect(const string &from, const string &to, int weight) {

  // if trying to make a self-edge, kill it
  if (from == to) {
    return false;
  }
  // check if both vertices exist
  if (!contains(from)) {
    add(from);
  }
  if (!contains(to)) {
    add(to);
  }
  // call helper corresponding to graph type
  if (directionalEdges) {
    return directionalConnection(from, to, weight);
  } else {
    return nondirectionalConnection(from, to, weight);
  }
}

/**
 * there can only be one edge between two vertices in a directed graph
 * so search the "from" vertex's adjVector for the "to" vertex then erase
 * @param from: the starting vertex
 * @param to: the destination vertex
 * @return bool whether or not the edge was removed successfully
 */
bool Graph::directionalDisconnection(const string &from, const string &to) {
  int eraseIndex = 0;
  // traverse the "from" vector until matching "to" value is found
  for (auto pair : theGraph.at(from)) {
    // when a match is found, erase that element & return true
    if (pair[1] == to) {
      theGraph[from].erase(theGraph.at(from).begin() + eraseIndex);
      return true;
    }
    eraseIndex++;
  }
  // if a match is never found return false
  return false;
}

/**
 * edges in a nondirectional graph must go both ways, so traverse both the
 * "from" and "to" vertices' adjVector and remove corresponding element
 * @param from: the starting vertex
 * @param to: the destination vertex
 * @return bool whether or not the edge was removed successfully
 */
bool Graph::nondirectionalDisconnect(const string &from, const string &to) {
  // variables tracking erase status & curr index
  bool firstErase = false;
  bool secondErase = false;
  int eraseIndex = 0;

  // erase connection from "from" vector
  for (auto pair : theGraph.at(from)) {
    if (pair[1] == to) {
      theGraph[from].erase(theGraph.at(from).begin() + eraseIndex);
      firstErase = true;
      break;
    }
    eraseIndex++;
  }

  // reset the index counter
  // erase connection from the "to" vector
  eraseIndex = 0;
  for (auto pair : theGraph.at(to)) {
    if (pair[1] == from) {
      theGraph[to].erase(theGraph.at(to).begin() + eraseIndex);
      secondErase = true;
      break;
    }
    eraseIndex++;
  }

  // return the results
  return (firstErase && secondErase);
}

/**
 * if either vertex does not exist return false, else call
 * the corresponding helper function for (non)directional graph
 * @param from: the starting vertex
 * @param to: the destination vertex
 * @return bool whether or not the edge was disconnected successfully
 */
bool Graph::disconnect(const string &from, const string &to) {
  // check if both vertices exist
  if ((!contains(from)) || (!contains(to))) {
    return false;
  }
  // disconnect based on graph type
  if (directionalEdges) {
    return directionalDisconnection(from, to);
  } else {
    return nondirectionalDisconnect(from, to);
  }
}

/**
 * recursive function which visits a node, prints it & adds to the results
 * string then checks if there is an unvisited vertex in the current vertex's
 * adjVector. if there is, recursive call on that vertex
 * @param results: the ongoing string containing results
 * @param visited: a map <string, bool> of visited nodes
 * @param startLabel: the vertex we are currently at
 */
string Graph::dfsHelper(string results, map<string, bool> &visited,
                        const string &startLabel) {
  visited[startLabel] = true;
  cout << startLabel;
  results += startLabel;
  vector<vector<string>> adjList = theGraph.at(startLabel);

  // should traverse in alphabetical order so
  // swap to do [0] label then [1] weight
  for (int i = 0; i < adjList.size(); i++) {
      string hold = adjList[i][0];
      adjList[i][0] = adjList[i][1];
      adjList[i][1] = hold;
  }
  // sort by the labels
  sort (adjList.begin(), adjList.end());
  // swap back to original [0] weight thn [1] label
  // implementing as a "vector of pairs" instead of a "vector of
  // vectors" would be better, be we're too far in now to go back >:)
  for (int i = 0; i < adjList.size(); i++) {
      string hold = adjList[i][0];
      adjList[i][0] = adjList[i][1];
      adjList[i][1] = hold;
  }

  // for each element in the adjacency list
  for (auto pair : adjList) {
    // if the number hasn't been visited
    if (visited.find(pair[1]) == visited.end()) {
      string newStart = pair[1];
      results = dfsHelper(results, visited, newStart);
    }
  }
  return results;
}

/**
 *
 *  *** DOES NOT PASS ASSERTS BUT DISPLAYS CORRECT ANSWERS ***
 *
 * if the starting vertex doesn't exist just return, else
 * call the dfs helper function to start recursive calls
 * @param startLabel: the starting vertex
 * @param visit: method to output traversal
 */
void Graph::dfs(const string &startLabel, void visit(const string &label)) {
  if (!contains(startLabel)) {
    return;
  }
  string results = "";
  map<string, bool> visited;
  results = dfsHelper(results, visited, startLabel);
  visit(results);
}

/**
 * Function that prints each node in a queue & adds it to result, then adds its
 * unvisited adjacent nodes to the queue (the breath-first approach to
 * traversing a graph).
 * @param results the string containing results.
 * @param visited a map <string, bool> of visited nodes.
 * @param startLabel the starting vertex of traversal.
 */
string Graph::bfsHelper(string results, map<string, bool> &visited,
                        const string &startLabel) {
  queue<string> q;
  q.push(startLabel);
  visited[startLabel] = true;

  while (!q.empty()) {
    string currentLabel = q.front();
    q.pop();
    cout << currentLabel;
    results += currentLabel;

    vector<vector<string>> adjList = theGraph.at(currentLabel);

    /* for each item in adjacency list, push unvisited nodes to queue and mark
    as visited */
    for (auto pair : adjList) {
      if (visited.find(pair[1]) == visited.end()) {
        q.push(pair[1]);
        visited[pair[1]] = true;
      }
    }
  }
  return results;
}

/**
 *
 *  *** DOES NOT PASS ASSERTS BUT DISPLAYS CORRECT ANSWERS ***
 *
 * if the starting vertex doesn't exist just return, else
 * call the bfs helper function to traverse graph
 * @param startLabel: the starting vertex
 * @param visit: method to output traversal
 */
void Graph::bfs(const string &startLabel, void visit(const string &label)) {
  if (!contains(startLabel)) {
    return;
  }

  string results = "";
  map<string, bool> visited;
  results = bfsHelper(results, visited, startLabel);

  visit(results);
}

/**
 * Finds shortest path between starting vertex and all other verticies in graph
 * using the Dijkstra algorithm approach. Stores the shortest distances in
 * map<string, int> weights, and the previous vertices in map<string, string>
 * previous.
 * @param startLabel the starting vertex to find shortest path from
 */
pair<map<string, int>, map<string, string>>
Graph::dijkstra(const string &startLabel) const {
  map<string, int> weights;
  map<string, string> previous;

  // priority queue to hold unvisited vertices
  typedef pair<int, string> p;
  priority_queue<p, vector<p>, greater<p>> pq;

  // if startLabel doesn't exist in graph, return
  if (!contains(startLabel)) {
    return make_pair(weights, previous);
  }

  // insert the starting key-val pair in graph and pq
  weights.insert({startLabel, 0});
  pq.push(make_pair(0, startLabel));

  while (!pq.empty()) {
    // get the top vertex, pop it
    string current = pq.top().second;
    pq.pop();

    vector<vector<string>> adjList = theGraph.at(current);

    // for all adjacent vertices of current vertex
    for (auto pair : adjList) {

      string vertex = pair[1];
      int weight = stoi(pair[0]);

      // calculate distance total distance to start vertex
      int dist = weight + weights[current];

      // if vertex not in weights/previous
      if (weights.find(vertex) == weights.end()) {
        weights.insert({vertex, dist});
        previous.insert({vertex, current});
        pq.push(make_pair(weights[vertex], vertex));
      }
      // if dist is less than current known dist of vertex
      else if (weights[vertex] > dist) {
        weights[vertex] = dist;
        previous[vertex] = current;
        pq.push(make_pair(weights[vertex], vertex));
      }
    }
  }

  weights.erase(startLabel); // erase the start vertex

  return make_pair(weights, previous);
}

// minimum spanning tree using Prim's algorithm
int Graph::mstPrim(const string &startLabel,
                   void visit(const string &from, const string &to,
                              int weight)) const {
  return -1;
}

// minimum spanning tree using Prim's algorithm
int Graph::mstKruskal(const string &startLabel,
                      void visit(const string &from, const string &to,
                                 int weight)) const {
  return -1;
}

// read a text file and create the graph
bool Graph::readFile(const string &filename) {
  ifstream myFile;
  myFile.open(filename);
  if (!myFile.is_open()) {
    cerr << "Failed to open " << filename << endl;
    return false;
  }
  int edges = 0;
  int weight = 0;
  string fromVertex;
  string toVertex;
  myFile >> edges;
  for (int i = 0; i < edges; ++i) {
    myFile >> fromVertex >> toVertex >> weight;
    connect(fromVertex, toVertex, weight);
  }
  myFile.close();
  return true;
}
