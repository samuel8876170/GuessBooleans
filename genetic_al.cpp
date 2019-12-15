#include "../include/genetic_al.h"

#define population_size 6
#define chromosome_size 10

int main(){
    Init();

    Answer ans(chromosome_size);
    ans.print();

    Population population(population_size, chromosome_size);
    population.print();

    int gen = 0;
    while (!population.IsCorrect(ans)) {
        printf("\n\nGeneration: %d\n", gen);

        population.roulette_wheel_selection(ans.getAnswer());

        population.cross_over();

        population.mutation(0.7);

        ++gen;
    }

    Chromosome result = population.getBestChromosome();
    result.print();

    return 0;
}
