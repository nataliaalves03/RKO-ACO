#ifndef _ACO_H
#define _ACO_H

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/************************************************************************************
 Method: UpdateArchiveSize
 Description: update the size of the archive solutions according the current state
*************************************************************************************/
static void UpdateArchiveSize(std::vector <TSol> &T, const int archive_size, const TProblemData &data)
{
    // size of the current archive
    int oldPsize = T.size();
    double best = T[0].ofv;

    // pruning 
    if (oldPsize > archive_size){
        T.resize(archive_size);
    }

    // generate new solutions 
    else if (oldPsize < archive_size){
        T.resize(archive_size);

        // Create the initial solutions with random keys 
        for (int i=oldPsize; i<archive_size; i++)
        {
            CreateInitialSolutions(T[i], data.n);
            T[i].ofv = Decoder(T[i], data);
        }

        std::sort(T.begin(), T.end(), sortByFitness);
    }
}


/************************************************************************************
 * Calcula o peso (omega_l) para uma solução no arquivo com base em seu rank. 
 * Utiliza a fórmula da Eq. (7) do artigo.
 *************************************************************************************/
double ACO_calculate_weight(const int archive_size, const double q_param, int rank) {
    // Implementa a Eq. (7) do artigo: omega_l = (1 / (qk*sqrt(2*pi))) * exp(-((l-1)^2) / (2 * q^2 * k^2))
    // O 'l' no artigo corresponde ao 'rank' aqui (1-based).
    double l_minus_1 = static_cast<double>(rank - 1);
    double qk_product = q_param * archive_size;
    
    // Evitar divisão por zero se qk_product for zero ou muito pequeno
    if (qk_product < std::numeric_limits<double>::epsilon()) {
        // Se qk_product é ~0, então o exp(-inf) é 0 exceto para l_minus_1 = 0.
        // Isso favorece fortemente o rank 1, como sugerido para 'q' pequeno.
        return (l_minus_1 == 0) ? 1.0 : 0.0;
    }

    double denominator = qk_product * std::sqrt(2.0 * M_PI);
    double exponent = -(l_minus_1 * l_minus_1) / (2.0 * qk_product * qk_product);
    return (1.0 / denominator) * std::exp(exponent);
}

/************************************************************************************
 * Calcula o desvio padrão (sigma_l^i) para uma dimensão e solução específica. 
 * Utiliza a fórmula da Eq. (9) do artigo.
 *************************************************************************************/
double ACO_calculate_std_dev(const int archive_size, const double xi_param, const std::vector<TSol>& T, int dimension_idx, int chosen_solution_idx) {
    // Implementa a Eq. (9) do artigo: sigma_l^i = xi * sum(|s_e^i - s_l^i|) / (k-1)
    // Onde 's_l' é a solução escolhida, 's_e' são as outras soluções no arquivo.
    double sum_distances = 0.0;
    int k_minus_1 = archive_size - 1;

    if (k_minus_1 == 0) { // Caso k=1, não há outras soluções para comparar
        return 0.0;
    }

    double s_l_i = T[chosen_solution_idx].rk[dimension_idx];

    for (int e = 0; e < archive_size; ++e) {
        if (e == chosen_solution_idx) continue; // Ignora a própria solução escolhida
        double s_e_i = T[e].rk[dimension_idx];
        sum_distances += std::abs(s_e_i - s_l_i);
    }

    return xi_param * (sum_distances / k_minus_1);
}

 /************************************************************************************
 * Constrói uma nova solução por uma "formiga". 
 * Envolve a amostragem de funções Gaussianas baseadas no arquivo.
 *************************************************************************************/
TSol ACO_construct_solution(const int archive_size, const double q_param, const double xi_param, const std::vector<TSol> &T, const int n) {
    TSol new_sol;
    new_sol.rk.resize(n);

    // 1. Escolher uma função Gaussiana (solução do arquivo) baseada nos pesos (Eq. 8)
    std::vector<double> weights(archive_size);
    double total_weight = 0.0;
    for (int l = 0; l < archive_size; ++l) {
        // rank é l+1 porque ranks são 1-baseados
        weights[l] = ACO_calculate_weight(archive_size, q_param, l + 1);
        total_weight += weights[l];
    }

    std::vector<double> probabilities(archive_size);
    for (int l = 0; l < archive_size; ++l) {
        probabilities[l] = weights[l] / total_weight;
    }

    // Amostragem ponderada para escolher uma solução do arquivo
    int chosen_solution_idx = -1;
    double rand_val = randomico (0, 1); 
    double cumulative_prob = 0.0;
    for (int l = 0; l < archive_size; ++l) {
        cumulative_prob += probabilities[l];
        if (rand_val <= cumulative_prob) {
            chosen_solution_idx = l;
            break;
        }
    }
    // Caso de borda, se rand_val for 1.0 ou por erro de ponto flutuante, escolher o último.
    if (chosen_solution_idx == -1) {
        chosen_solution_idx = archive_size - 1;
    }

    // 2. Amostrar a função Gaussiana escolhida para cada random key
    for (int d = 0; d < n; ++d) {
        double mu_val = T[chosen_solution_idx].rk[d];
        double sigma_val = ACO_calculate_std_dev(archive_size, xi_param, T, d, chosen_solution_idx);

        // Garantir que sigma_val não seja zero ou muito pequeno para evitar problemas na distribuição normal
        if (sigma_val < std::numeric_limits<double>::epsilon()) {
            sigma_val = 0.9999; // Usar uma dispersão padrão grande
        }

        std::normal_distribution<double> normal_dist(mu_val, sigma_val);
        double sampled_value = normal_dist(rng);

        // rejection sampling
        while (sampled_value < 0 || sampled_value >= 1)
        {
            sampled_value = normal_dist(rng);
        }

        new_sol.rk[d] = sampled_value;
    }

    return new_sol;
}

/************************************************************************************
 * Atualiza o arquivo de soluções com novas soluções geradas. 
 * Mantém apenas as k melhores soluções.
 *************************************************************************************/
void ACO_update_archive(const int archive_size, std::vector<TSol> &T, const std::vector<TSol>& new_solutions) {
    // Adiciona as novas soluções ao arquivo temporariamente
    for (const auto& sol : new_solutions) {
        T.push_back(sol);
    }

    // Re-ordena todo o arquivo
    std::sort(T.begin(), T.end(), sortByFitness);

    // Remove as piores soluções até que o arquivo tenha o tamanho archive_size
    if ((int)T.size() > archive_size) {
        T.resize(archive_size);
    }
}

/************************************************************************************
 Method: ACO
 Description: search process of the Ant Colony Optimization
*************************************************************************************/
void ACO(const TRunData &runData, const TProblemData &data)
{
    const char* method = "ACO";

    int archive_size = 0;                        // tamanho do arquivo de soluções (k).
    int num_ants     = 0;                        // número de formigas (soluções) geradas por iteração (m).
    double q_param   = 0;                        // parâmetro 'q' que influencia a diversificação/intensificação.
    double xi_param  = 0;                        // parâmetro 'xi' que influencia a velocidade de convergência.

    std::vector<TSol> T;                         // O arquivo de soluções (memória de feromônio)
    TSol bestSolution;                           // A melhor solução global encontrada até agora

    // double bestOFcurrent = 0;                // best ofv found in the current generation

    // local variables
    int numGenerations = 0;                  // number of generations
    float currentTime = 0;                   // computational time of the search process
    int bestGeneration = 0;                  // number of generation that found the best solution
    int improv = 0;                          // improvement flag

    double start_timeMH = get_time_in_seconds();    // start computational time
    double end_timeMH = get_time_in_seconds();      // end computational time

    std::vector<int> RKorder;                   // define a order for the neighors
    RKorder.resize(data.n);
    std::iota(RKorder.begin(), RKorder.end(), 0);

    // Q-Learning parameters
    std::vector<TState> S;                      // finite state space
    int numPar = 0;                             // number of parameters
    int numStates = 0;                          // number of states
    int iCurr = 0;                              // current (initial) state
    double epsilon=0;                           // greed choice possibility
    double lf=0;                                // learning factor
    double df=0;                                // discount factor
    double R=0;                                 // reward
    std::vector <std::vector <TQ> > Q;          // Q-Table
    std::vector<int> ai;                        // actions
    float epsilon_max = 1.0;                    // maximum epsilon 
    float epsilon_min = 0.1;                    // minimum epsilon
    int Ti = 1;                                 // number of epochs performed
    int restartEpsilon = 1;                     // number of restart epsilon
    int st = 0;                                 // current state
    int at = 0;                                 // current action

    // ** read file with parameter values
    numPar = 4;
    std::vector<std::vector<double>> parameters;
    parameters.resize(numPar);

    readParameters(method, runData.control, parameters, numPar);

    // offline control
    if (runData.control == 0){
        // define parameters of ACO (offline)
        archive_size = parameters[0][0];
        num_ants     = parameters[1][0];
        q_param      = parameters[2][0];
        xi_param     = parameters[3][0];
    }

    // online control
    else{
        // Q-Learning 
        if (runData.control == 1){
            // create possible states of the Markov chain
            CreateStates(parameters, numStates, numPar, S);

            // number of restart epsilon
            restartEpsilon = 1;  

            // maximum epsilon  
            epsilon_max = 1.0;  

            // current state
            iCurr = irandomico(0,numStates-1);

            // define parameters of ACO
            archive_size = (int)S[iCurr].par[0];
            num_ants = S[iCurr].par[1];
            q_param = S[iCurr].par[2];
            xi_param  = S[iCurr].par[3];
        }
    }  

    // initialize solution archive T
    T.clear();
    T.resize(archive_size);

    bestSolution.ofv = INFINITY;
    for (int i = 0; i < archive_size; ++i) {
        // initialize solutions
        
        CreateSolutionSemiGreedy(T[i], data);
        
        if (T[i].ofv < bestSolution.ofv) {
            bestSolution = T[i];
        }
    }
    // Sort archive T in increase order of fitness
    std::sort(T.begin(), T.end(), sortByFitness);

    UpdatePoolSolutions(bestSolution, method, runData.debug);

    double bestSolutionGen;
    
    // run the evolutionary process until stop criterion
    while (currentTime < runData.MAXTIME*runData.restart && !stop_execution.load())
    {
    	// number of generations
        numGenerations++;
        bestSolutionGen = INFINITY;

        // define the parameters considering the current state and evolve a new iteration of the ACO
        // Q-Learning 
        if (runData.control == 1){
            // set Q-Learning parameters  
            SetQLParameter(currentTime, Ti, restartEpsilon, epsilon_max, epsilon_min, epsilon, lf, df, runData.MAXTIME*runData.restart); 

            // choose a action at for current state st
            at = ChooseAction(S, st, epsilon);

            // execute action at of st
            iCurr = S[st].Ai[at];

            // define the parameters according of the current state
            num_ants     = S[iCurr].par[1];
            q_param      = S[iCurr].par[2];
            xi_param     = S[iCurr].par[3];

            if(archive_size != (int)S[iCurr].par[0]){
                archive_size = (int)S[iCurr].par[0];
                UpdateArchiveSize(T, archive_size, data);
            }
        }

        std::vector<TSol> newSolutions;
        newSolutions.reserve(num_ants);

        // Fase de Construção de Soluções por Formigas
        for (int i = 0; i < num_ants; ++i) 
        {    
            if (stop_execution.load()) return;   
            
            TSol newSol = ACO_construct_solution(archive_size, q_param, xi_param, T, data.n);
            newSol.ofv = Decoder(newSol, data); 
            newSolutions.push_back(newSol);

            // Atualiza a melhor solução 
            if (newSol.ofv < bestSolution.ofv) {
                bestSolution = newSol;
                bestGeneration = numGenerations; 
                improv = 1;
                UpdatePoolSolutions(bestSolution, method, runData.debug);
            }
            if(newSol.ofv < bestSolutionGen){
                bestSolutionGen = newSol.ofv;
            }
        }

        // Fase de Atualização do Feromônio (Arquivo de Soluções)
        ACO_update_archive(archive_size, T, newSolutions);

        // local search
        TSol newSol = bestSolution;
        NelderMeadSearch(newSol, data);
        
        // Atualiza a melhor solução 
        if (newSol.ofv < bestSolution.ofv) {
            bestSolution = newSol;
            bestGeneration = numGenerations; 
            improv = 1;
            UpdatePoolSolutions(bestSolution, method, runData.debug);
        }
        if(newSol.ofv < bestSolutionGen){
            bestSolutionGen = newSol.ofv;
        }

        //if (runData.debug) printf("\nGen: %d [%d, %d, %.2lf, %.2lf] \t sBest: %lf ", numGenerations, archive_size, num_ants, q_param, xi_param, bestSolution.ofv);

        // Q-Learning 
        if (runData.control == 1){
            // The reward function is based on improvement of the current best fitness and binary reward
            if (improv){
                R = 1;
                improv = 0;
            }
            else{
                R = (bestSolution.ofv - bestSolutionGen)/bestSolutionGen;
            }

            // if (runData.debug) printf("\t [%.4lf, %d, %d] \t [%d]", R, st, at, Psize);

            // index of the next state
            int st_1 = S[st].Ai[at];

            // Update the Q-Table value
            S[st].Qa[at] = S[st].Qa[at] + lf*(R + df*S[st_1].maxQ - S[st].Qa[at]); 

            if (S[st].Qa[at] > S[st].maxQ)
            {
                S[st].maxQ = S[st].Qa[at];
                S[st].maxA = at;
            }

            // Define the new current state st
            st = st_1;
        }
        
        // terminate the evolutionary process in MAXTIME
        end_timeMH = get_time_in_seconds();
        currentTime = end_timeMH - start_timeMH;
    }

    // free memory of ACO components
    T.clear();
    bestSolution.rk.clear();

    // print policy
    // if (runData.debug and runData.control == 1)
    //     PrintPolicy(S, st);
}

#endif