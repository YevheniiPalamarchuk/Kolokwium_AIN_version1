using System;
using System.Collections.Generic;
using System.Linq;

namespace Kolokwium_AIN_version1
{
    internal class Program
    {
        // Individual
        public class Individual
        {
            public double[] Genes { get; set; }
            public double Fitness { get; set; }
        }

        // List 
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

        // Evaluate population
        public static void EvaluatePopulation(List<Individual> population, Func<double[], double> fitnessFunction)
        {
            foreach (var individual in population)
            {
                individual.Fitness = fitnessFunction(individual.Genes);
            }
        }
        // GeneralizedRosenbrock
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

        //CalculateFitnessGeneralizedRosenbrock
        public static double CalculateFitnessGeneralizedRosenbrock(double[] genes)
        {
            double a = 1.0;
            double b = 100.0;

            double result = GeneralizedRosenbrockFunction(genes, a, b);

            // Since minimizing, negate the result
            double fitness = -result;

            return fitness;
        }

        // Salomon
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

        // Whitley
        public static double CalculateFitnessWhitley(double[] genes)
        {
            int n = genes.Length;
            double sum = 0;

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    double term1 = 100 * Math.Pow(genes[i] * genes[i] - genes[j], 2);
                    double term2 = Math.Pow(1 - genes[j], 2);
                    sum += (term1 + term2) / 4000 - Math.Cos(term1 + term2) + 1;
                }
            }

            return sum;
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


        public static List<Individual> UniformCrossover(List<Individual> parents, double crossoverRate)
        {
            List<Individual> offspring = new List<Individual>();

            Random rand = new Random();

            for (int i = 0; i < parents.Count; i += 2)
            {
                if (rand.NextDouble() < crossoverRate)
                {
                    double[] child1Genes = new double[parents[i].Genes.Length];
                    double[] child2Genes = new double[parents[i].Genes.Length];

                    for (int j = 0; j < parents[i].Genes.Length && j < parents[i + 1].Genes.Length; j++)
                    {
                        // Perform crossover for valid indices
                        if (rand.NextDouble() < 0.5)
                        {
                            child1Genes[j] = parents[i].Genes[j];
                            child2Genes[j] = parents[i + 1].Genes[j];
                        }
                        else
                        {
                            child1Genes[j] = parents[i + 1].Genes[j];
                            child2Genes[j] = parents[i].Genes[j];
                        }
                    }

                    offspring.Add(new Individual { Genes = child1Genes });
                    offspring.Add(new Individual { Genes = child2Genes });
                }
                else
                {
                    // If no crossover, simply copy parents to offspring
                    offspring.Add(new Individual { Genes = parents[i].Genes.ToArray() });
                    offspring.Add(new Individual { Genes = parents[i + 1].Genes.ToArray() });
                }
            }

            return offspring;
        }



        public static void Mutate(List<Individual> population, double mutationRate, double mutationRange)
        {
            Random rnd = new Random();

            foreach (var individual in population)
            {
                if (rnd.NextDouble() < mutationRate)
                {
                    for (int i = 0; i < individual.Genes.Length; i++)
                    {
                        // Add a small random value to the gene
                        individual.Genes[i] += rnd.NextDouble() * 2 * mutationRange - mutationRange;
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

            int populationSize = 100;
            int numVariables = 3;
            double minValueGeneralizedRosenbrock = -30.0;
            double maxValueGeneralizedRosenbrock = 30.0;

            double minValueSalomon = -100.0;
            double maxValueSalomon = 100.0;

            double minValueWhitley = -10.24;
            double maxValueWhitley = 10.24;

            int numGenerations = 100;

            // Dimensions
            int[] dimensions = { 10, 30, 50 };

            //List<Individual> population = InitializePopulation(populationSize, numVariables, minValue, maxValue);
            List<Individual> populationGeneralizedRosenbrock = InitializePopulation(populationSize, numVariables, minValueGeneralizedRosenbrock, maxValueGeneralizedRosenbrock);


            for (int generation = 0; generation < numGenerations; generation++)
            {
                // Evaluate the fitness of each individual
                EvaluatePopulation(populationGeneralizedRosenbrock, CalculateFitnessGeneralizedRosenbrock); // or CalculateFitnessGeneralizedRosenbrock

                // Select parents for the next generation
                int numParents = populationSize / 2;
                List<Individual> selectedParents = RouletteWheelSelection(populationGeneralizedRosenbrock, numParents);

                // Perform crossover on selected parents
                double crossoverRate = 0.7;
                List<Individual> offspring = UniformCrossover(selectedParents, crossoverRate);

                // Perform mutation on the offspring
                double mutationRate = 0.1;
                double mutationRange = 0.5;
                Mutate(offspring, mutationRate, mutationRange);

                // Replace the old population with the new one using generational replacement
                ReplaceGenerational(populationGeneralizedRosenbrock, offspring);

                // Optionally, print or log information about the current generation
                Console.WriteLine($"Generation {generation + 1}: Best Fitness = {populationGeneralizedRosenbrock.Min(ind => ind.Fitness)}");
            }

            List<Individual> populationSalomon = InitializePopulation(populationSize, numVariables, minValueSalomon, maxValueSalomon);

            for (int generation = 0; generation < numGenerations; generation++)
            {
                // Evaluate the fitness of each individual using the Salomon function
                EvaluatePopulation(populationSalomon, CalculateFitnessSalomon);

                // Select parents for the next generation
                int numParents = populationSize / 2; // Adjust as needed
                List<Individual> selectedParents = RouletteWheelSelection(populationSalomon, numParents);

                // Perform crossover on selected parents
                double crossoverRate = 0.7;
                List<Individual> offspring = UniformCrossover(selectedParents, crossoverRate);

                // Perform mutation on the offspring
                double mutationRate = 0.1;
                double mutationRange = 0.5;
                Mutate(offspring, mutationRate, mutationRange);

                // Replace the old population with the new one using generational replacement
                ReplaceGenerational(populationSalomon, offspring);

                // Optionally, print or log information about the current generation
                Console.WriteLine($"Generation {generation + 1} (Salomon): Best Fitness = {populationSalomon.Min(ind => ind.Fitness)}");
            }

            List<Individual> populationWhitley = InitializePopulation(populationSize, numVariables, minValueWhitley, maxValueWhitley);

            for (int generation = 0; generation < numGenerations; generation++)
            {
                // Evaluate the fitness of each individual using the Whitley function
                EvaluatePopulation(populationWhitley, CalculateFitnessWhitley);

                // Select parents for the next generation
                int numParents = populationSize / 2; // Adjust as needed
                List<Individual> selectedParents = RouletteWheelSelection(populationWhitley, numParents);

                // Perform crossover on selected parents
                double crossoverRate = 0.7;
                List<Individual> offspring = UniformCrossover(selectedParents, crossoverRate);

                // Perform mutation on the offspring
                double mutationRate = 0.1;
                double mutationRange = 0.5;
                Mutate(offspring, mutationRate, mutationRange);

                // Replace the old population with the new one using generational replacement
                ReplaceGenerational(populationWhitley, offspring);

                // Optionally, print or log information about the current generation
                Console.WriteLine($"Generation {generation + 1} (Whitley): Best Fitness = {populationWhitley.Min(ind => ind.Fitness)}");
            }
        }
    }
}