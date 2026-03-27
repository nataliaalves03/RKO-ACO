// *******************************************************************
//      file with specific functions to solve the problem
// *******************************************************************
#ifndef _PROBLEM_QMC_VSBPP_H
#define _PROBLEM_QMC_VSBPP_H

#include <unordered_map>
#include <unordered_set>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <sstream>
#include <algorithm>
#include <random>


#include <list>
#include <mutex>

// Cache for decoded vectors and their FO values
class FO_Cache {
public:
    using Key = std::string;
    using Value = double;
    static const size_t MAX_SIZE = 1000;

    // Insert or update cache
    void put(const Key& key, Value value) {
        std::lock_guard<std::mutex> lock(mtx);

        auto it = map.find(key);
        if (it != map.end()) {
            // Move to back (most recent)
            order.erase(it->second.second);
            order.push_back(key);
            it->second = {value, --order.end()};
        } else {
            // New entry
            if (order.size() >= MAX_SIZE) {
                // Remove oldest
                const Key& oldest = order.front();
                map.erase(oldest);
                order.pop_front();
            }
            order.push_back(key);
            map[key] = {value, --order.end()};
        }
    }

    // Check if key exists
    bool get(const Key& key, Value& value) {
        std::lock_guard<std::mutex> lock(mtx);

        auto it = map.find(key);
        if (it != map.end()) {
            value = it->second.first;
            return true;
        }
        return false;
    }

private:
    std::list<Key> order;
    std::unordered_map<Key, std::pair<Value, std::list<Key>::iterator>> map;
    std::mutex mtx;
};

static FO_Cache fo_cache;



//----------------- DEFINITION OF PROBLEM SPECIFIC TYPES -----------------------


struct Bin {
    int cost;
    std::vector<int> capacities;
    int total_capacity;
};

struct Item {
    std::vector<int> weights;
    int total_link_cost;
    int total_weight;
    int total_link_weight;
    std::vector<int> links;
};

struct TProblemData
{
    int n;      // size of the RKO vector 

    // other variables of the problem at hand
    int N; // number of items
    int M; // number of bins
    int D; // number of dimensions
    std::vector<Item> items;
    std::vector<Bin> bins;
    std::vector<std::vector<int>> links;
    std::vector<std::vector<int>> graph;
    std::vector<int> binsByCapacity;
    std::vector<int> itemsByLinkCost;
    std::vector<int> itemsByLinkWeight;
    int qtd_params;
    int default_bin_type;
};



//----------------------- IMPLEMENTATION OF FUNCTIONS  -------------------------------

void ReadData(char name[], TProblemData &data);
void PreProcess_graph(TProblemData &data);
void CreateSolutionSemiGreedy(TSol &s, const TProblemData &data);

double calculateObjective(const TProblemData& data, const std::unordered_map<int, UsedBin>& usedBins);
UsedBin openNewBin(const TProblemData& data, int bin_type);
bool placeItemToBin(const TProblemData& data, int item_id, UsedBin &usedBin);
bool AssignItemToBestBin(const TProblemData& data, std::unordered_map<int, UsedBin>& usedBins, int& binCounter, int i);
void AssignItemsToBins(const TSol& s, const std::vector<int>& sol, const TProblemData& data, std::unordered_map<int, UsedBin>& usedBins);
void AssignSemiGreedy(const std::vector<int>& sol, const TProblemData& data, std::unordered_map<int, UsedBin>& usedBins, std::vector<bool>& assigned, int RCL_size, int& binCounter);
void reduceBinCost(const TProblemData& data, UsedBin& bin);
void reduceUsedBinsCost(const TProblemData& data, std::unordered_map<int, UsedBin>& usedBins);
void mergeBins(const TProblemData& data, std::unordered_map<int, UsedBin>& usedBins);
void itemRelocationLS(const TProblemData& data, std::unordered_map<int, UsedBin>& usedBins, double p, int max_iter);
void localSearch(const TProblemData& data, std::unordered_map<int, UsedBin>& usedBins, TSol &s);

void DecodeRk(TSol &s, const TProblemData &data, std::vector<int>& sol);
void EncodeRk(TSol &s, const TProblemData &data, const std::unordered_map<int, UsedBin>& usedBins);
double Decoder(TSol &s, const TProblemData &data);

void FreeMemoryProblem(TProblemData &data);
void PrintSolution(TSol& s, const TProblemData& data, const char *filename);
void printBins(TSol& s, const TProblemData &data, double fo, std::unordered_map<int, UsedBin> usedBins);


/************************************************************************************
 Method: ReadData
 Description: read input data of the problem
*************************************************************************************/
void ReadData(char name[], TProblemData &data)
{ 
    FILE *arq;
    arq = fopen(name,"r");

    if (arq == NULL){
        printf("\nERROR: File (%s) not found!\n",name);
        getchar();
        exit(1);
    }

    char line[1024];
    bool found_instance = false;

    // 1. Find the start of the next instance
    while (fgets(line, sizeof(line), arq))
    {
        if (strncmp(line, "instance=", 9) == 0)
        {
            found_instance = true;
            break;
        }
    }
    if (!found_instance)
    {
        fclose(arq);
        printf("No instance found in file.\n");
        exit(1);
    }

    // 2. Parse instance header
    sscanf(line, "instance=%*d,n=%d,m=%d,d=%d", &data.N, &data.M, &data.D);

    data.items.clear();
    data.bins.clear();
    data.items.reserve(data.N);
    data.bins.reserve(data.M);
    data.links.clear();
    data.graph.clear();
    // N x N matrix
    data.links.resize(data.N);
    data.graph.resize(data.N);
    for (int i = 0; i < data.N; ++i) {
        data.links[i].resize(data.N, 0);
        data.graph[i].resize(data.N, 0);
    }


    // 3. Parse bin types (sorted by cost)
    for (int i = 0; i < data.M; ++i) {
        // Read BinType line
        if (!fgets(line, sizeof(line), arq))
            break;

        Bin bin;
        sscanf(line, "BinType:%*d, cost=%d", &bin.cost);
        
        // Read d lines for bin capacities
        int cap = 0;
        for (int j = 0; j < data.D; ++j)
        {
            if (!fgets(line, sizeof(line), arq))
                break;
            
            sscanf(line, "    no=%*d, value=%d", &cap);
            bin.capacities.push_back(cap);
        }
        data.bins.push_back(bin);
    }

    // 4. Parse items
    for (int i = 0; i < data.N; ++i) {
        // Read Item line
        if (!fgets(line, sizeof(line), arq))
            break;
        // Read d lines for item weights
        Item item;
        int weight = 0;
        for (int j = 0; j < data.D; ++j)
        {
            if (!fgets(line, sizeof(line), arq))
                break;
           
            sscanf(line, "   no=%*d, value=%d", &weight);
            item.weights.push_back(weight);
        }
        data.items.push_back(item);
    }

    // 5. Skip to the end of the instance (line with slashes)
    while (fgets(line, sizeof(line), arq)) {
        if (strstr(line, "////") != NULL) 
            break;
    }

    // 6. Check for "Links between items" and read links
    if (fgets(line, sizeof(line), arq) && strstr(line, "Links between items") != NULL) {
        // read all link lines
        while (fgets(line, sizeof(line), arq)) {

            int from = 0, to = 0, cost = 0;
            if (sscanf(line, "Link %*d: itemNo1=%d, itemNo2=%d, cost=%d", &from, &to, &cost) == 3) {
                //zero index
                from = from-1;
                to = to-1;
                data.links[from][to] = cost;
                data.links[to][from] = cost;
            }
        }
    }

    fclose(arq);



    // ---- Pre processamento

    PreProcess_graph(data);


    //get order of bins by aggragated capacity
    data.binsByCapacity.resize(data.M, 0);
    
    for (int m = 0; m < data.M; m++) {
        data.binsByCapacity[m] = m;
        data.bins[m].total_capacity = std::accumulate(data.bins[m].capacities.begin(), data.bins[m].capacities.end(), 0);
    }
    std::sort(data.binsByCapacity.begin(), data.binsByCapacity.end(), [&data](int a, int b) {
        return data.bins[a].total_capacity < data.bins[b].total_capacity;
    });

    data.default_bin_type = data.binsByCapacity[data.M-1]; 


    //get items total weight
    for (int i = 0; i < data.N; i++) {
        data.items[i].total_weight = std::accumulate(data.items[i].weights.begin(), data.items[i].weights.end(), 0);
    }

    //get items links weight
    for (int i = 0; i < data.N; i++) {
        int link_weight = 0;
        for(int j = 0; j < data.N; j++){
            if(data.links[i][j] > 0){
                link_weight += data.items[j].total_weight;
            }
        }
        data.items[i].total_link_weight = link_weight;
    }

    data.itemsByLinkCost.resize(data.N, 0);
    data.itemsByLinkWeight.resize(data.N, 0);

    for(int i = 0; i < data.N; i++) {
        data.itemsByLinkCost[i] = i;
        data.itemsByLinkWeight[i] = i;
    }
    //descending
    std::sort(data.itemsByLinkCost.begin(), data.itemsByLinkCost.end(), [&data](int a, int b) {
        return data.items[a].total_link_cost > data.items[b].total_link_cost;
    });
    std::sort(data.itemsByLinkWeight.begin(), data.itemsByLinkWeight.end(), [&data](int a, int b) {
        return data.items[a].total_link_weight > data.items[b].total_link_weight;
    });

    //List of links
    for (int i = 0; i < data.N; i++) {
        Item item = data.items[i];
        for (int j = 0; j < data.N; j++) {
            if(data.links[i][j] > 0 && data.graph[i][j] != -1){
                item.links.push_back(j);
            }
        }
    }

    data.qtd_params = 3;
    data.n = data.N + data.qtd_params; //Params: max bins, alloc strategy
   
    //exit(1);
}


void PreProcess_graph(TProblemData &data){

    data.graph = data.links;
   
    for (int i = 0; i < data.N; ++i) {
        for (int j = i+1; j < data.N; ++j) {
            if(data.graph[i][j] >= 0){
                bool fit = true;
                
                for(int m = 0; m < data.M; m++){
                    auto& item_i = data.items[i];
                    auto& item_j = data.items[j];

                    for(int d = 0; d < data.D; d++){
                        if(item_i.weights[d] + item_j.weights[d] > data.bins[m].capacities[d]){
                            fit = false;
                            break;
                        }
                    }
                    if(fit) {
                        item_i.total_link_cost += data.links[i][j];
                        item_j.total_link_cost += data.links[j][i];
                        break;
                    }
                }
                //só há link se os items podem ser adicionados juntos em um bin.
                if (!fit) {
                    data.graph[i][j] = -1;
                    data.graph[j][i] = -1;
                    //printf("link quebrado: %d-%d\n",i,j);
                }
            }
        }
    }

}




// -----------------------
// Create Initial Solutions
// -----------------------

void CreateSolutionSemiGreedy(TSol &s, const TProblemData &data) {
    s.rk.resize(data.n);

    std::vector<bool> assigned(data.N, false);
    std::unordered_map<int, UsedBin> usedBins;
    int binCounter = 0;
    
    //RCL size is between 3 and 5% of N
    int RCL_size = std::uniform_int_distribution<int>(3, (int)(data.N * 0.05))(rng);

    std::vector<int> sol(data.n, 0);
    for (int i = 0; i < data.N; i++) {
        sol[i] = i; // allocation order 0..N-1
    }

    AssignSemiGreedy(sol, data, usedBins, assigned, RCL_size, binCounter);

    //assign the rest of items
    for(int i = 0; i < data.N; i++){
        if(assigned[i]) continue;
        assigned[i] = AssignItemToBestBin(data, usedBins, binCounter, i);
    }
    
    reduceUsedBinsCost(data, usedBins);
    mergeBins(data, usedBins);
    itemRelocationLS(data, usedBins, 1, 500);
    mergeBins(data, usedBins);
    EncodeRk(s, data, usedBins);
    s.ofv = calculateObjective(data, usedBins);
    s.usedBins = usedBins;
}



// -----------------------
// Post process and Local Search
// -----------------------

void mergeBins(const TProblemData& data, std::unordered_map<int, UsedBin>& usedBins) {
    // Track which bins have been merged to avoid overlap
    std::unordered_set<int> mergedBins;
    std::vector<std::tuple<int, int, UsedBin>> merges; // (binIdx_a, binIdx_b, newBin)

    for (auto it_a = usedBins.begin(); it_a != usedBins.end(); ++it_a) {
        int binIdx_a = it_a->first;
        UsedBin& bin_a = it_a->second;
        if (bin_a.items.empty() || mergedBins.count(binIdx_a)) continue;

        for (auto it_b = usedBins.begin(); it_b != usedBins.end(); ++it_b) {
            int binIdx_b = it_b->first;
            UsedBin& bin_b = it_b->second;
            if (bin_b.items.empty() || binIdx_a >= binIdx_b || mergedBins.count(binIdx_b)) continue;

            int bins_cost = data.bins[bin_a.binType].cost + data.bins[bin_b.binType].cost;
            if(data.bins[data.M-1].cost > bins_cost) break;

            for (int m = 0; m < data.M; m++) {
                Bin m_bin = data.bins[m];
                if(m_bin.cost > bins_cost) break;
                bool fit = true;
                for (int d = 0; d < data.D; d++) {
                    if (bin_a.usedCapacity[d] + bin_b.usedCapacity[d] > m_bin.capacities[d]) {
                        fit = false;
                        break;
                    }
                }
                if (fit) {
                    UsedBin newBin = openNewBin(data, m);
                    for (int d = 0; d < data.D; d++) {
                        newBin.usedCapacity[d] = bin_a.usedCapacity[d] + bin_b.usedCapacity[d];
                    }
                    newBin.items.insert(newBin.items.end(), bin_a.items.begin(), bin_a.items.end());
                    newBin.items.insert(newBin.items.end(), bin_b.items.begin(), bin_b.items.end());
                    merges.emplace_back(binIdx_a, binIdx_b, newBin);
                    mergedBins.insert(binIdx_a);
                    mergedBins.insert(binIdx_b);
                    usedBins[binIdx_a].items.clear();
                    usedBins[binIdx_b].items.clear();
                    break; // Only merge each pair once
                }
            }
            if (mergedBins.count(binIdx_a)) break; // bin_a already merged
        }
    }

    // Remove merged bins and add new bins
    int newIdx = 0;
    for (const auto& pair : usedBins) {
        if (pair.first >= newIdx) newIdx = pair.first + 1;
    }

    for (const auto& [a, b, newBin] : merges) {
        //usedBins.erase(a);
        //usedBins.erase(b);
        usedBins[newIdx++] = newBin;
    }
}



void reduceBinCost(const TProblemData& data, UsedBin& bin) {
    for (int m = 0; m < bin.binType; m++) { //bins sorted by cost
        bool fits = true;
        for (int d = 0; d < data.D; d++) {
            if (data.bins[m].capacities[d] <= bin.usedCapacity[d]) {
                fits = false;
                break;
            }
        }
        if (fits) {
            bin.binType = m;
            break;
        }
    }
}

void reduceUsedBinsCost(const TProblemData& data, std::unordered_map<int, UsedBin>& usedBins) {
    for (auto& [binIdx, bin] : usedBins) {
        if (bin.items.empty()) continue;
        reduceBinCost(data, bin);
    }
}



// Move items to reduce link/bin cost and possibly close bins
void itemRelocationLS(const TProblemData& data, std::unordered_map<int, UsedBin>& usedBins, double p, int max_iter) {
    bool improved = true;
    //double fo = calculateObjective(data, usedBins);
    
    struct ItemMoveCandidate {
        int item;
        int fromBin;
        int externalLinkCost;
        int binCost;
        int binSize;
    };

    while (improved && max_iter > 0) {
        improved = false;

        // 1. Compute priority for each item: link cost to items in other bins
        std::vector<ItemMoveCandidate> candidates;
        for (const auto& [binIdx, bin] : usedBins) {
            if (bin.items.empty()) continue;
            double r = std::uniform_real_distribution<double>(0, 1)(rng);
            if(r <= p){
                for (int item : bin.items) {
                    candidates.push_back({
                        item,
                        binIdx,
                        data.items[item].total_link_weight,
                        data.bins[bin.binType].cost,
                        (int)bin.items.size()
                    });
                }
            }
        }

        // Sort: prioritize high external link cost, then expensive bins, then small bins
        std::sort(candidates.begin(), candidates.end(), [](const ItemMoveCandidate& a, const ItemMoveCandidate& b) {
            if (a.externalLinkCost != b.externalLinkCost)
                return a.externalLinkCost > b.externalLinkCost;
            if (a.binCost != b.binCost)
                return a.binCost > b.binCost;
            return a.binSize < b.binSize;
        });

        // 2. Try to move each candidate item to another bin
        for (const auto& cand : candidates) {
            int item = cand.item;
            int fromBin = cand.fromBin;

            // Try all other bins
            for (auto& [toBinIdx, toBin] : usedBins) {
                if (toBinIdx == fromBin) continue;
                if (max_iter-- <= 0) break;

                // Check if item fits in toBin
                UsedBin tempBin = toBin;
                if (!placeItemToBin(data, item, tempBin)) continue;

                // Simulate move: remove from fromBin, add to toBin
                UsedBin tempFromBin = usedBins[fromBin];
                // Remove item from tempFromBin
                auto it = std::find(tempFromBin.items.begin(), tempFromBin.items.end(), item);
                if (it != tempFromBin.items.end()) {
                    for (int d = 0; d < data.D; ++d)
                        tempFromBin.usedCapacity[d] -= data.items[item].weights[d];
                    tempFromBin.items.erase(it);
                }

                reduceBinCost(data, tempBin);
                if(!tempFromBin.items.empty()) reduceBinCost(data, tempFromBin);

                //Evaluate new solution
                double deltaObjective = 0.0;

                // 1. Bin cost delta
                double oldFromBinCost = data.bins[usedBins[fromBin].binType].cost;
                double oldToBinCost   = data.bins[usedBins[toBinIdx].binType].cost;
                double newFromBinCost = tempFromBin.items.empty() ? 0 : data.bins[tempFromBin.binType].cost;
                double newToBinCost   = data.bins[tempBin.binType].cost;
                deltaObjective += (newFromBinCost + newToBinCost) - (oldFromBinCost + oldToBinCost);

                // 2. Link cost delta
                // For all items in fromBin and toBin, recalculate only the part affected by the moved item
                for (int other : usedBins[fromBin].items) {
                    if (other == item) continue;
                    deltaObjective += data.links[item][other]; // They were together, now apart
                }
                for (int other : usedBins[toBinIdx].items) {
                    deltaObjective -= data.links[item][other]; // They were apart, now together
                }

                if(deltaObjective < 0){
                    usedBins[toBinIdx] = std::move(tempBin);
                    usedBins[fromBin] = std::move(tempFromBin);
                    improved = true;
                    break;
                }
                
            }
            if (improved || max_iter <= 0) break;
        }
    }
}



void localSearch(const TProblemData& data, std::unordered_map<int, UsedBin>& usedBins, TSol &s){
    double p = s.rk[data.n-3];
    if(p >= 0.8){
        itemRelocationLS(data, usedBins, p, 100);
        mergeBins(data, usedBins);
    }
}

// -----------------------
// Calculate Objective Function
// -----------------------

double calculateObjective(const TProblemData& data, const std::unordered_map<int, UsedBin>& usedBins){
    bool feasible = true;
    double fo = 0.0;
    std::vector<int> itemsAllocated(data.N, 0);
    std::vector<int> itemsToBins(data.N, -1);

    for (auto& [binIdx, bin] : usedBins) {
        if(bin.items.empty()) continue;

        const Bin& bin_data = data.bins[bin.binType];
        int bin_cost = bin_data.cost;

        //bin cost
        fo += bin_cost;

        //capacities
        for(int d = 0; d < data.D; d++){
            int bin_cap = bin_data.capacities[d];
            int used_cap = bin.usedCapacity[d];

            if(bin_cap < used_cap){
                //penalty infeasible by capacity
                fo += used_cap - bin_cap;
                feasible = false;
                printf("infeasible by capacity\n");
                exit(1);
            }
            else{
                //desempate. mais vazio, melhor
                fo += (double) used_cap / bin_cap / 1000000.0;

                //desempate. melhor custo baneficio
                fo += (((double) bin_cost / bin_cap) + (1.0 / bin_cap)) / 100000.0;
            }
        }

        //items
        for(int i : bin.items){
            itemsAllocated[i]++;
            itemsToBins[i] = binIdx;
        }
    }


    //links cost
    double link_cost = 0;
    double link_solved = 0;
    for(int i = 0; i < data.N; i++){
        for(int j = i+1; j < data.N; j++){
            if(data.links[i][j] > 0 && itemsToBins[i] != -1 && itemsToBins[j] != -1 && itemsToBins[i] != itemsToBins[j]){
                link_cost += data.links[i][j];
            }
            if(data.links[i][j] > 0 && itemsToBins[i] != -1 && itemsToBins[j] != -1 && itemsToBins[i] == itemsToBins[j]){
                //reward solved links
                link_solved += data.links[i][j];
            }
        }
    }
    fo += link_cost;
    fo += link_cost / (link_solved + 1e-8) / 1e4;


    return fo;
}



// -----------------------
// Assign items to bins and set bin type
// -----------------------


UsedBin openNewBin(const TProblemData& data, int bin_type){
    UsedBin newBin;
    newBin.binType = bin_type;
    newBin.items.clear();
    newBin.usedCapacity.resize(data.D, 0);
    return newBin;
}


bool placeItemToBin(const TProblemData& data, int item_id, UsedBin& usedBin){
    Item item = data.items[item_id];
    Bin bin = data.bins[usedBin.binType];

    //Sum check
    int total_used_capacity = 0;
    for(int d = 0; d < data.D; d++){
        total_used_capacity +=  usedBin.usedCapacity[d];
    }
    if(item.total_weight + total_used_capacity > bin.total_capacity){
        return false;
    }
    
    //Check feasibility per d
    for(int d = 0; d < data.D; d++){
        if(item.weights[d] + usedBin.usedCapacity[d] > bin.capacities[d])
            return false; //exceeded capacity in d.
    }

    //Place item
    for(int d = 0; d < data.D; d++){
        usedBin.usedCapacity[d] += item.weights[d];
    }
    usedBin.items.push_back(item_id);

    return true;
}

void AssignSemiGreedy(const std::vector<int>& sol, const TProblemData& data, std::unordered_map<int, UsedBin>& usedBins, std::vector<bool>& assigned, int RCL_size, int& binCounter) {
    int unassignedCount = data.N;
    while (unassignedCount > 0) {
        int start = -1;

        // 1. Create a RCL
        std::vector<int> RCL;
        for(int j = 0; j < data.N && RCL.size() < RCL_size; j++){
            //allocation from the rk/sol order
            int i = sol[j];
            if(!assigned[i]) RCL.push_back(i);
        }
        if (RCL.empty()) break;
        // Get the unassigned item with higher total_link_weight
        std::sort(RCL.begin(), RCL.end(), [&data](int a, int b) {
            return data.items[a].total_link_weight > data.items[b].total_link_weight;
        });
        start = RCL.front();
        
        if (start == -1) break;

        // 2. Place the start item
        UsedBin bin = openNewBin(data, data.default_bin_type);
        placeItemToBin(data, start, bin);
        assigned[start] = true;
        unassignedCount--;
        
        // 3. Greedily add closest unassigned items by link cost
        while (true) {
            int best = -1, best_cost = 0;
            for (int i = 0; i < data.N; ++i) {
                if (!assigned[i] && data.links[start][i] > 0 && data.graph[start][i] != -1 && data.links[start][i] > best_cost) {
                    UsedBin tempBin = bin;
                    if (placeItemToBin(data, i, tempBin)) {
                        best = i;
                        best_cost = data.links[start][i];
                    }
                }
            }
            if (best == -1) break;
            placeItemToBin(data, best, bin);
            assigned[best] = true;
            unassignedCount--;
        }
        usedBins[binCounter++] = bin;
    }
}

bool AssignItemToBestBin(const TProblemData& data, std::unordered_map<int, UsedBin>& usedBins, int& binCounter, int i) {
    Item item_i = data.items[i];
    int link_cost = item_i.total_link_cost;

    int bestBin = -1;
    int bestLinkCost = link_cost;
    bool placed = false;
    int new_link_cost = 0;
 
    //find best existing bin (link cost) that fits
    for (const auto& [binIdx, bin] : usedBins) {
        
        //emulate link cost
        new_link_cost = link_cost;
        for(int bin_item : bin.items){
            new_link_cost -= data.links[i][bin_item];
        }
        
        //check capacity
        if(new_link_cost < bestLinkCost){
            UsedBin temp = bin;
            if(placeItemToBin(data, i, temp)){
                bestBin = binIdx;
                bestLinkCost = new_link_cost;
            }
        }
    }

    //if found, add item to best bin
    if( bestBin != -1){
        placed = placeItemToBin(data, i, usedBins[bestBin]);
    }

    //otherwise, open new bin
    if(!placed){
        UsedBin new_bin = openNewBin(data, data.default_bin_type);
        if(placeItemToBin(data, i, new_bin)){
            placed = true;
            usedBins[binCounter++] = new_bin;
        }
    }

    return placed;
}

void AssignItemsToBins(const TSol& s, const std::vector<int>& sol, const TProblemData& data, std::unordered_map<int, UsedBin>& usedBins) {
    usedBins.clear();
    int binCounter = 0;
    std::vector<bool> assigned(data.N, false);
    int param_max_bins = round(data.N * s.rk[data.n-1]);
    int param_strategy = sol[data.n-2];

    if(param_strategy == 0){
        //semi-greedy
        int RCL_size = std::max(3, (int)(data.N * 0.05));
        AssignSemiGreedy(sol, data, usedBins, assigned, RCL_size, binCounter);
    }

    //Param = 1: open max bins with 1 item each / Param 2: skip it.
    if(param_strategy == 1){
        for(int b = 0; b < param_max_bins; b++){
            int i = sol[b];
            UsedBin new_bin = openNewBin(data, data.default_bin_type);
            if(placeItemToBin(data, i, new_bin)){
                assigned[i] = true;
                usedBins[binCounter++] = new_bin;
            }
        }
    }

    if(param_strategy == 1 || param_strategy == 2){
            
        //assign the items with link cost
        for(int s = 0; s < data.N; s++){
            int i = sol[s];
            if(assigned[i]) continue;

            Item item_i = data.items[i];
            int link_cost = item_i.total_link_cost;

            if(link_cost == 0) continue; //no link cost, assign later
            
            assigned[i] = AssignItemToBestBin(data, usedBins, binCounter, i);
        }

    }
    

    //assign the rest of items
    for(int s = 0; s < data.N; s++){
        int i = sol[s];
        if(assigned[i]) continue;
        assigned[i] = AssignItemToBestBin(data, usedBins, binCounter, i);
    }
}





/************************************************************************************
 Method: Decoders
 Description: mapping the random-key solutions into problem solutions
*************************************************************************************/

void DecodeRk(TSol &s, const TProblemData &data, std::vector<int>& sol){ 
    for (int j = 0; j < data.N; j++) sol[j] = j;
    std::sort(sol.begin(), sol.begin() + data.N, [&s](int i1, int i2) {
        return s.rk[i1] < s.rk[i2];
    });

    //Param max bins
    sol[data.n-1] = round(data.N * s.rk[data.n-1]);

    //Param assign strategy
    int qtd_strategies = 3;
    sol[data.n-2] = std::min(static_cast<int>(s.rk[data.n-2] * qtd_strategies + 1e-8), qtd_strategies - 1);

    //Param item relocation p
    sol[data.n-3] = s.rk[data.n-3] >= 0.7;
}



// Updates s.rk based on the final usedBins assignment
void EncodeRk(TSol &s, const TProblemData &data, const std::unordered_map<int, UsedBin>& usedBins) {
    std::vector<int> itemOrder(data.N, -1);
    int idx = 0;

    for (const auto& [binIdx, bin] : usedBins) {
        for (int item : bin.items) {
            itemOrder[idx++] = item;
        }
    }

    for (int i = 0; i < data.N; ++i) {
        s.rk[i] = (double)itemOrder[i] / data.N - 0.000001; // [0, 1) 
    }

    //Atualiza params
    int countBins = 0;
    for(auto [idx,bin] : usedBins){
        if(!bin.items.empty()) countBins++;
    }
    s.rk[data.n-1] = countBins / data.N;

    //Param assign strategy (initial sol. settings)
    s.rk[data.n-2] = 0.9; // strategy=2

    //Param item relocation p
    s.rk[data.n-3] = 0.8; // enabled

}

double Decoder(TSol &s, const TProblemData &data) {  
    std::vector<int> sol(data.n);
    DecodeRk(s, data, sol);

    // Serialize sol to string as cache key
    std::ostringstream oss;
    for (int v : sol) oss << v << ",";
    std::string key = oss.str();

    double cached_fo;
    if (fo_cache.get(key, cached_fo)) {
        s.ofv = cached_fo;
        return s.ofv;
    }

    std::unordered_map<int, UsedBin> usedBins;
    AssignItemsToBins(s, sol, data, usedBins);

    reduceUsedBinsCost(data, usedBins);
    mergeBins(data, usedBins);

    localSearch(data, usedBins, s);

    s.ofv = calculateObjective(data, usedBins);
    s.usedBins = usedBins;

    // Store in cache
    fo_cache.put(key, s.ofv);
    
    return s.ofv;
}




/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem(TProblemData &data)
{
    //
}


void printBins(TSol& s, const TProblemData &data, double fo, std::unordered_map<int, UsedBin> usedBins){
         
    printf("\n--- Solution (used bins) ---\n");

    bool feasible = true;

    // Check feasibility of all used bins
    for (const auto& [binIdx, bin] : usedBins) {
        if (bin.items.empty()) continue;
        const Bin& bin_data = data.bins[bin.binType];

        // Check capacity constraints for each dimension
        for (int d = 0; d < data.D; ++d) {
            if (bin.usedCapacity[d] > bin_data.capacities[d]) {
                feasible = false;
                break;
            }
        }
        if (!feasible) break;
    }

    if (!feasible) {
        printf("WARNING: Solution is infeasible (bin capacity exceeded).\n");
    }


    int bins_cost = 0;

    for (const auto& [binIdx, bin] : usedBins) {
        if (bin.items.empty()) continue;
        const Bin& bin_data = data.bins[bin.binType];

        printf("Bin %2d (type %d): Used [", binIdx, bin.binType);
        for (int d = 0; d < data.D; ++d) {
            printf("%ld", bin.usedCapacity[d]);
            if (d < data.D - 1) printf(", ");
        }
        printf("] / [");
        for (int d = 0; d < data.D; ++d) {
            printf("%d", bin_data.capacities[d]);
            if (d < data.D - 1) printf(", ");
        }
        printf("], Cost %d | Items:", bin_data.cost);
        bins_cost +=  bin_data.cost;

        for (int item_id : bin.items) {
            const Item& item = data.items[item_id];
            printf(" %d(", item_id);
            for (int d = 0; d < data.D; ++d) {
                printf("%d", item.weights[d]);
                if (d < data.D - 1) printf(",");
            }
            printf(")");
        }
        printf("\n");
    }
    printf("\n");


    std::vector<int> itemAssignmentCount(data.N, 0);
    for (const auto& [binIdx, bin] : usedBins) {
        for (int item : bin.items) {
            itemAssignmentCount[item]++;
        }
    }
    for (int i = 0; i < data.N; ++i) {
        if (itemAssignmentCount[i] == 0) {
            printf( "WARNING: Item %d is not assigned to any bin!\n", i);
        }
        if (itemAssignmentCount[i] > 1) {
            printf( "WARNING: Item %d is assigned to multiple bins (%d times)!\n", i, itemAssignmentCount[i]);
        }
    }


    for (const auto& [binIdx, bin] : usedBins) {
        if (bin.items.empty()) continue;
        const Bin& bin_data = data.bins[bin.binType];
        for (int d = 0; d < data.D; ++d) {
            if (bin.usedCapacity[d] > bin_data.capacities[d]) {
                printf( "WARNING: Bin %d (type %d) exceeds capacity in dimension %d: used %ld > cap %d\n",
                    binIdx, bin.binType, d, bin.usedCapacity[d], bin_data.capacities[d]);
            }
        }
    }


    // Print links cost
    std::vector<int> itemToBin(data.N, -1);
    for (const auto& [binIdx, bin] : usedBins) {
        for (int item : bin.items) {
            itemToBin[item] = binIdx;
        }
    }
    int links_cost = 0;
    for (int i = 0; i < data.N; ++i) {
        for (int j = i + 1; j < data.N; ++j) {
            if (data.links[i][j] > 0 && itemToBin[i] != -1 && itemToBin[j] != -1 && itemToBin[i] != itemToBin[j]) {
                links_cost += data.links[i][j];
            }
        }
    }
    printf("Links cost: %d\n", links_cost);
    printf("Bins cost: %d\n", bins_cost);
    fo =  links_cost + bins_cost;

    printf("Calculated FO: %4.f \t Solution FO: %4.f  \n", fo, s.ofv);
    
}

void PrintSolution(TSol& s, const TProblemData& data, const char *filename) {

    std::vector<int> sol(data.n);
    DecodeRk(s, data, sol);
    std::unordered_map<int, UsedBin> usedBins;
    usedBins = s.usedBins;
    double fo = calculateObjective(data, usedBins);
    

    //print to file
    if (strcmp(filename,"") != 0) {

        FILE *solFile;
        solFile = fopen(filename,"a");

        fprintf(solFile,"\n--- Solution (used bins) ---\n");

        bool feasible = true;

        // Check feasibility of all used bins
        for (const auto& [binIdx, bin] : usedBins) {
            if (bin.items.empty()) continue;
            const Bin& bin_data = data.bins[bin.binType];

            // Check capacity constraints for each dimension
            for (int d = 0; d < data.D; ++d) {
                if (bin.usedCapacity[d] > bin_data.capacities[d]) {
                    feasible = false;
                    break;
                }
            }
            if (!feasible) break;
        }

        if (!feasible) {
            fprintf(solFile,"WARNING: Solution is infeasible (bin capacity exceeded).\n");
        }


        int bins_cost = 0;

        for (const auto& [binIdx, bin] : usedBins) {
            if (bin.items.empty()) continue;
            const Bin& bin_data = data.bins[bin.binType];

            fprintf(solFile,"Bin %2d (type %d): Used [", binIdx, bin.binType);
            for (int d = 0; d < data.D; ++d) {
                fprintf(solFile,"%ld", bin.usedCapacity[d]);
                if (d < data.D - 1) fprintf(solFile,", ");
            }
            fprintf(solFile,"] / [");
            for (int d = 0; d < data.D; ++d) {
                fprintf(solFile,"%d", bin_data.capacities[d]);
                if (d < data.D - 1) fprintf(solFile,", ");
            }
            fprintf(solFile,"], Cost %d | Items:", bin_data.cost);
            bins_cost +=  bin_data.cost;

            for (int item_id : bin.items) {
                const Item& item = data.items[item_id];
                fprintf(solFile," %d(", item_id);
                for (int d = 0; d < data.D; ++d) {
                    fprintf(solFile,"%d", item.weights[d]);
                    if (d < data.D - 1) fprintf(solFile,",");
                }
                fprintf(solFile,")");
            }
            fprintf(solFile,"\n");
        }
        fprintf(solFile,"\n");


        // Print links cost
        std::vector<int> itemToBin(data.N, -1);
        for (const auto& [binIdx, bin] : usedBins) {
            for (int item : bin.items) {
                itemToBin[item] = binIdx;
            }
        }
        int links_cost = 0;
        for (int i = 0; i < data.N; ++i) {
            for (int j = i + 1; j < data.N; ++j) {
                if (data.links[i][j] > 0 && itemToBin[i] != -1 && itemToBin[j] != -1 && itemToBin[i] != itemToBin[j]) {
                    links_cost += data.links[i][j];
                }
            }
        }
        fprintf(solFile,"Links cost: %d\n", links_cost);
        fprintf(solFile,"Bins cost: %d\n", bins_cost);


        fprintf(solFile,"Calculated FO: %4.f \t Solution FO: %4.f  \n", fo, s.ofv);
       

        fclose(solFile);
    }
    else { 

       for (int i = 0; i < sol.size(); ++i) {
           printf("%d ", sol[i]);
       }
       printf("\n");

       //printBins(s, data, fo, usedBins);
    }
}


#endif
