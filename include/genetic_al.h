#ifndef C__PROJECT_GENETIC_AL_H
#define C__PROJECT_GENETIC_AL_H

#endif //C__PROJECT_GENETIC_AL_H
#include <iostream>
//#include <bits/stdc++.h>
#include <vector>
#include <random>
#include <ctime>

using namespace std;

void Init() {
    // initialize random seed first
    srand((unsigned int) time(nullptr));
}


class Answer{
private:
    vector<bool> answer;
public:
    // randomly generate N boolean
    Answer(int N){
        answer.clear();
        for (int i = 0; i < N; ++i)
            answer.push_back(rand()%10 > 5);
    }



//////////////////////////////output function/////////////////////////////////
    vector<bool>& getAnswer(){
        return answer;
    }
    const void print(){
        printf("answer: [");
        for (unsigned int i = 0; i < answer.size()-1; ++i)
            printf("%s, ", answer[i]?"1":"0");
        printf("%s]\n", *(answer.end()-1)?"1":"0");
    }
};


class Chromosome{
private:
    vector<bool> chromosome;
    int chromosome_size;
public:
    // randomly generate N gene
    Chromosome(int N){
        chromosome.clear();
        chromosome_size = N;
        for (int i = 0; i < N; ++i)
            chromosome.push_back(rand()%10 > 5?1:0);
    }

    int fitness(vector<bool> &answer){
        if (chromosome_size != answer.size())
            printf("Error(fitness): chromosome.size != answer.size!\n");

        int fitness = 0;
        for (int i = 0; i < chromosome_size; ++i) {
            if (chromosome[i] == answer[i])
                ++fitness;
        }
        return fitness;
    }

    void cross_over(Chromosome &other){
        // Step 2: randomly choose random number of genes to crossover from first chromosome
        // get a number of genes which need crossover
        int cross_gen_num = 1 + rand() % chromosome_size;

        vector<int> cross_index;
        int i = 0, tempi;
        // get a list of genes which need crossover
//        printf("cross_index: ");
        gen_index : for (; i < cross_gen_num; ++i) {
            tempi = rand() % chromosome_size;
            for (unsigned int j = 0; j < cross_index.size(); ++j) {
                if (tempi == cross_index[j])
                    goto gen_index;
            }
            cross_index.push_back(tempi);
//            printf("%d, ", tempi);
        }
//        printf("\n");

        sort(cross_index.begin(), cross_index.end());   // sort before swap genes for convenience
        bool tempb;
        for (int j = 0, k = 0; j < chromosome_size && k < cross_index.size(); ++j) {
            if (j == cross_index[k]) {
                tempb = chromosome[j];
                chromosome[j] = other.chromosome[j];
                other.chromosome[j] = tempb;
                ++k;
            }
        }
    }

    void mutate(){
        // Step 2: randomly choose "10%" of genes to mutate (at least 2)
        // Step 3: change their order randomly
        int genes_num = 2 +  (int)(0.1 * chromosome_size), i = 0, tempi;

        vector<int> gene_index, backup;
        gene_index.reserve(genes_num);
        backup.reserve(genes_num);

        // Step 2
        gen_index : for (; i < genes_num; ++i) {
            tempi = rand() % chromosome_size;
            for (unsigned int j = 0; j < gene_index.size(); ++j) {
                if (gene_index[j] == tempi)
                    goto gen_index;
            }
            gene_index.push_back(tempi);
        }
        backup = gene_index;


        // Step 3
        i = 0;
        while(gene_index.size() > 0){
            tempi = rand() % gene_index.size();
            chromosome[backup[i]] = gene_index[tempi];
            gene_index.erase(gene_index.begin() + tempi);
            ++i;
        }
    }



//////////////////////////////output function/////////////////////////////////
    vector<bool>& getChromosome(){
        return chromosome;
    }
    const void print(){
        printf("[");
        for (int i = 0; i < chromosome_size-1; ++i)
            printf("%s, ", chromosome[i]?"1":"0");
        printf("%s]\n", *(chromosome.end()-1)?"1":"0");
    }
};


class Population{
private:
    vector<Chromosome> population;
    vector<double> cumulative_fitness_prob;
    int population_size;
    int most_fit_index;

public:
    // generate populationSize chromosome
    Population(int populationSize, int chromosome_size) {
        if (populationSize < 2)
            printf("The population size less than 2 which will have error!!\n");

        population_size = populationSize;
        cumulative_fitness_prob.reserve(populationSize);

        population.clear();
        for (int i = 0; i < populationSize; ++i)
            population.push_back(Chromosome(chromosome_size));
    }

    void fitness_probability(vector<bool> &answer){
        cumulative_fitness_prob.clear();
        int sum = 0, temp, min = 2147483647, max = -2147483647;

        vector<double> fitness_prob;
//        printf("fitness : ");
        for (int i = 0; i < population_size; ++i) {
            temp = population[i].fitness(answer);

            fitness_prob.push_back(temp);
            sum += temp;
//            printf("%d, ", temp);
        }

        // print fitness prob
//        printf("fitness_prob : [\n");
        for (int i = 0; i < population_size; ++i) {
            if (fitness_prob[i] > 0)
                fitness_prob[i] /= sum;
//            printf("%.3f, ", fitness_prob[i]);
        }
//        printf("]\n\n");

        // make it cumulative below~
        cumulative_fitness_prob.push_back(0);
        double s = 0;
        for (int i = 0; i < population_size; ++i) {
            s += fitness_prob[i];
            cumulative_fitness_prob.push_back(s);
        }

        // print cumulative fitness prob
//        printf("cumulative_fitness_prob : [\n");
//        for (unsigned int i = 0; i < cumulative_fitness_prob.size()-1; ++i)
//            printf("%.3f, ", cumulative_fitness_prob[i]);
//        printf("%.3f]\n\n", *(cumulative_fitness_prob.end()-1));
    }

    void roulette_wheel_selection(vector<bool> &answer){
        fitness_probability(answer);

        vector<Chromosome> new_population;
        vector<double> prob;
        for (int i = 0; i < population_size; ++i)
            prob.push_back((double)(rand() % RAND_MAX) / RAND_MAX);

        sort(prob.begin(), prob.end());

        for (int i = 0; i < population_size; ++i) {
            for (unsigned int j = 1; j < cumulative_fitness_prob.size(); ++j) {
                if (prob[i] < cumulative_fitness_prob[j] && prob[i] >= cumulative_fitness_prob[j-1]) {
                    new_population.push_back(population[j - 1]);
                    break;
                }
            }
        }

        if (population_size != new_population.size())
            printf("Error: population.size():%d != new_population.size():%d\n",
                    population.size(), new_population.size());
        population = new_population;
    }

    void cross_over(){
        // Step 1: randomly choose random number of chromosome as parents (at least 2)
        // Step 2: then randomly choose random number of genes to crossover from first chromosome
        int parent_num = 2 + (rand() % (population_size-1)) / 2 * 2;
        vector<int> parent_index;
        int temp, i = 0;
        gen_index : for (; i < parent_num; ++i) {
            temp = rand() % population_size;
            for (unsigned int j = 0; j < parent_index.size(); ++j) {
                if (temp == parent_index[j])
                    goto gen_index;
            }
            parent_index.push_back(temp);        // Step 1
        }

        //// Debug
        if (parent_index.size() != parent_num)
            printf("Error: parent_index.size():%d != parent_num:%d\n", parent_index.size(), parent_num);
//        printf("Parent_index: ");
//        for (i = 0; i < parent_num; ++i){
//            printf("%d, ", parent_index[i]);
//        }
//        printf("\n");
        //// Debug

        for (i = 0; i < parent_num-1; i += 2)
            // Step 2
            population[parent_index[i]].cross_over(population[parent_index[i+1]]);
    }
    
    void mutation(double threshold_prob){
        // threshold_prob should [0, 1] and high enough

        // Step 1: randomly choose random number of chromosome to do mutation (with a threshold_prob test)
        // Step 2: randomly choose "10%" of genes to mutate (at least 2)
        // Step 3: change their order randomly
        int mutate_init_num = 1 + rand() % population_size, i = 0, temp;
        vector<int> mutate_index;
        gen_index : for (; i < mutate_init_num; ++i) {
            temp = rand() % population_size;
            for (unsigned int j = 0; j < mutate_index.size(); ++j) {
                if (temp == mutate_index[j])
                    goto gen_index;
            }
            mutate_index.push_back(temp);
        }

        //// Debug
        if (mutate_index.size() != mutate_init_num)
            printf("Error: mutate_index.size():%d != mutate_init_num:%d\n", mutate_index.size(), mutate_init_num);
//        printf("Parent_index: ");
//        for (i = 0; i < mutate_init_num; ++i){
//            printf("%d, ", mutate_index[i]);
//        }
//        printf("\n");
        //// Debug

        for (unsigned int j = 0; j < mutate_index.size(); ++j) {
            if ((double)rand() / RAND_MAX > threshold_prob)
                population[mutate_index[j]].mutate();
        }
    }

    bool IsCorrect(Answer &answer){
        vector<bool> ans = answer.getAnswer();
        printf("fitness: ");
        for (int i = 0; i < population_size; ++i) {
            if (population[i].fitness(ans) == ans.size()) {
                printf("\"%d\", ", population[i].fitness(ans));
                most_fit_index = i;
                return true;
            } else {
                printf("%d, ", population[i].fitness(ans));
            }
        }
        printf("\n");
        return false;
    }



//////////////////////////////output function/////////////////////////////////
    Chromosome getBestChromosome(){
        printf("most_fit_index : %d\n", most_fit_index);
        return population[most_fit_index];
    }
    const void print(){
        printf("population: [\n");
        for (int i = 0; i < population_size-1; ++i)
            population[i].print();
        (population.end()-1)->print();
        printf("]\n\n");
    }
};
