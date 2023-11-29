using System;
using System.Collections.Generic;

namespace Kolokwium_AIN_version1
{
    internal class Program
    {
        public class Individual
        {
            public double[] Genes { get; set; }
            public double Fitness { get; set; }
        }
        public List<Individual> InitializePopulation(int populationSize, int numVariables, double minValue, double maxValue)
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
        public void EvaluatePopulation(List<Individual> population)
        {
            foreach (var individual in population)
            {
                individual.Fitness = CalculateFitness(individual.Genes);
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
        public double CalculateFitness(double[] genes)
        {
            // Implement the fitness calculation for your specific function
            // For example, use the Generalized Rosenbrock function
            // double fitness = // ... calculate fitness based on the function ...
            double fitness = 0;
            return fitness;
        }

        static void Main(string[] args)
        {
            Console.WriteLine("Hello World!");
        }
    }
}
