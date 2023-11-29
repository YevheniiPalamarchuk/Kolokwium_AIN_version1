using System;
using System.Collections.Generic;
using System.Linq;

namespace Kolokwium_AIN_version1
{
    internal class Program
    {
        public class Individual
        {
            public double[] Genes { get; set; }
            public double Fitness { get; set; }
        }
        public static List<Individual> InitializePopulation(int populationSize, int numVariables, double minValue, double maxValue)
        {
            Random rand = new Random();
            List<Individual> population = new List<Individual>();

            for (int i = 0; i < populationSize; i++)
            {
                double[] genes = new double[numVariables];
                for (int j = 0; j < numVariables; j++)
                {
                    genes[j] = rand.NextDouble() * (maxValue - minValue) + minValue;
                }

                population.Add(new Individual { Genes = genes });
            }

            return population;
        }
        public static void EvaluatePopulation(List<Individual> population, Func<double[], double> fitnessFunction)
        {
            foreach (var individual in population)
            {
                individual.Fitness = fitnessFunction(individual.Genes);
            }
        }

        static double GeneralizedRosenbrockFunction(double[] x, double a, double b)
        {
            int n = x.Length;
            double sum = 0;

            for (int i = 0; i < n - 1; i++)
            {
                double term1 = a * Math.Pow(x[i + 1] - Math.Pow(x[i], 2), 2);
                double term2 = b * Math.Pow(1 - x[i], 2);
                sum += term1 + term2;
            }

            return sum;
        }
        public static double CalculateFitnessGeneralizedRosenbrock(double[] genes)
        {
            double a = 1.0;
            double b = 100.0;

            double result = GeneralizedRosenbrockFunction(genes, a, b);

            // Since you're minimizing, negate the result
            double fitness = -result;

            return fitness;
        }

        public static double CalculateFitnessSalomon(double[] genes)
        {
            double sum = 0;

            foreach (var gene in genes)
            {
                sum += Math.Pow(gene, 2);
            }

            double sqrtSum = Math.Sqrt(sum);

            return 1 - Math.Cos(2 * Math.PI * sqrtSum) + 0.1 * sqrtSum;
        }

        public static List<Individual> RouletteWheelSelection(List<Individual> population, int numParents)
        {
            List<Individual> parents = new List<Individual>();

            // Calculate total fitness of the population
            double totalFitness = population.Sum(individual => individual.Fitness);

            // Create a roulette wheel
            double[] wheel = new double[population.Count];
            double cumulativeFitness = 0;

            // Assign slots on the wheel based on individual fitness
            for (int i = 0; i < population.Count; i++)
            {
                wheel[i] = cumulativeFitness + population[i].Fitness / totalFitness;
                cumulativeFitness = wheel[i];
            }

            // Spin the wheel to select parents
            for (int i = 0; i < numParents; i++)
            {
                double spin = new Random().NextDouble();
                int selectedParentIndex = Array.FindIndex(wheel, value => value >= spin);

                parents.Add(population[selectedParentIndex]);
            }

            return parents;
        }


        public static List<Individual> OnePointCrossover(List<Individual> parents, double crossoverRate)
        {
            List<Individual> offspring = new List<Individual>();

            for (int i = 0; i < parents.Count; i += 2)
            {
                if (new Random().NextDouble() < crossoverRate)
                {
                    int crossoverPoint = new Random().Next(parents[i].Genes.Length);

                    double[] child1Genes = parents[i].Genes.Take(crossoverPoint)
                        .Concat(parents[i + 1].Genes.Skip(crossoverPoint))
                        .ToArray();

                    double[] child2Genes = parents[i + 1].Genes.Take(crossoverPoint)
                        .Concat(parents[i].Genes.Skip(crossoverPoint))
                        .ToArray();

                    offspring.Add(new Individual { Genes = child1Genes });
                    offspring.Add(new Individual { Genes = child2Genes });
                }
                else
                {
                    // If no crossover, offspring are identical to parents
                    offspring.Add(new Individual { Genes = parents[i].Genes });
                    offspring.Add(new Individual { Genes = parents[i + 1].Genes });
                }
            }

            return offspring;
        }
        public static void Mutate(List<Individual> population, double mutationRate, double mutationRange)
        {
            Random rand = new Random();

            foreach (var individual in population)
            {
                if (rand.NextDouble() < mutationRate)
                {
                    for (int i = 0; i < individual.Genes.Length; i++)
                    {
                        // Add a small random value to the gene
                        individual.Genes[i] += rand.NextDouble() * 2 * mutationRange - mutationRange;
                    }
                }
            }
        }

        public static void ReplaceGenerational(List<Individual> population, List<Individual> offspring)
        {
            // Combine the current population and the offspring
            List<Individual> combinedPopulation = new List<Individual>(population);
            combinedPopulation.AddRange(offspring);

            // Sort the combined population by fitness (assuming lower fitness is better)
            combinedPopulation.Sort((a, b) => a.Fitness.CompareTo(b.Fitness));

            // Select the top N individuals to form the new population
            for (int i = 0; i < population.Count; i++)
            {
                population[i] = combinedPopulation[i];
            }
        }


        // Body
        static void Main(string[] args)
        {

            int populationSize = 50;
            int numVariables = 3;
            double minValue = -30.0; 
            double maxValue = 30.0; 
            int numGenerations = 100;

            List<Individual> population = InitializePopulation(populationSize, numVariables, minValue, maxValue);

            for (int generation = 0; generation < numGenerations; generation++)
            {
                // Evaluate the fitness of each individual
                EvaluatePopulation(population, CalculateFitnessSalomon); // or CalculateFitnessGeneralizedRosenbrock

                // Select parents for the next generation
                int numParents = populationSize / 2; 
                List<Individual> selectedParents = RouletteWheelSelection(population, numParents);

                // Perform crossover on selected parents
                double crossoverRate = 0.7;
                List<Individual> offspring = OnePointCrossover(selectedParents, crossoverRate);

                // Perform mutation on the offspring
                double mutationRate = 0.1;
                double mutationRange = 0.5;
                Mutate(offspring, mutationRate, mutationRange);

                // Replace the old population with the new one using generational replacement
                ReplaceGenerational(population, offspring);

                // Optionally, print or log information about the current generation
                Console.WriteLine($"Generation {generation + 1}: Best Fitness = {population.Min(ind => ind.Fitness)}");
            }
        }
    }
}
