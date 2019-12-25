#pragma once
#include "Gen_Alg.h"

std::ofstream out{ "resultRK.txt" };

double fitness(std::pair<double, double> unit)
{
    double fit = 1 / (1 + (unit.first * unit.first) + (unit.second * unit.second));
    return fit;
}

std::vector<std::pair<double, double>> Generator()
{
    std::vector<std::pair<double, double>> unitVector(4);
    for (int i = 0; i < 4; i++)
    {
        std::mt19937 engine(std::random_device{}());
        auto generator = std::uniform_real_distribution<double>(0.0, 4.0);
        double x = generator(engine) - 2;
        generator = std::uniform_real_distribution<double>(0.0, 4.0);
        double y = generator(engine) - 2;
        unitVector[i] = std::make_pair(x, y);
    }

    std::vector<double> fit(4);
    out << "0(start)";
    double maxFit = -1;
    double sumFit = 0;
    for (int i = 0; i < 4; i++)
    {
        fit[i] = fitness(unitVector[i]);
        sumFit += fitness(unitVector[i]);
    }
    for (int i = 0; i < 4; i++)
    {
        maxFit = std::max(maxFit, fit[i]);
    }
    out << "\t|" << std::fixed << std::setprecision(5) << unitVector[0].first
        << " \t|" << unitVector[0].second << " \t|" << fit[0] << " \t|" << std::setw(16) << " |\n";
    for (int i = 1; i < 4; i++)
    {
        out << "\t\t|" << std::fixed << std::setprecision(5) << unitVector[i].first
            << " \t|" << unitVector[i].second << " \t|" << fit[i] << " \t|" << std::setw(16) << "|\n";
    }
    out << "\t\t|" << std::setw(16) << "|" << std::setw(16) << "|" << std::setw(16) << "|" << maxFit << std::setw(8) << "|" << sumFit / 4 << "\n";
    out << "---------------------------------------------------------------------------------------\n";
    return unitVector;
}

std::pair<std::pair<double, double>, std::pair<double, double>> children(std::pair<double, double> parent1,
                                                                         std::pair<double, double> parent2)
{
    std::pair<double, double> child1;
    child1.first = parent1.first;
    child1.second = parent2.second;
    std::pair<double, double> child2;
    child2.first = parent2.first;
    child2.second = parent1.second;
    return std::make_pair(child1, child2);
}

std::pair<double, double> mutation(std::pair<double, double> unit)
{
    std::mt19937 engine(std::random_device{}());
    int parameter = rand() % 2;
    auto generatorDouble = std::uniform_real_distribution<double>(0.0, 0.2);
    double delta = generatorDouble(engine) - 0.1;
    if (parameter == 0)
    {
        unit.first += delta;
    }
    else
    {
        unit.second += delta;
    }
    return unit;
}

std::pair<int, int> parentGenerator(std::vector<double> fitness, double sumFit)
{
    std::pair<int, int> parent;
    std::mt19937 engine(std::random_device{}());
    auto generator = std::uniform_real_distribution<double>(0.0, sumFit);
    double p1 = generator(engine);
    int chosenIndex;
    int index1 = 0, index2 = 0;
    bool foundCorrect = false;
    while (!foundCorrect)
    {
        double sum = 0;
        index1++;
        for (int j = 0; j < index1; j++)
        {
            sum += fitness[j];
            if (p1 <= sum)
            {
                foundCorrect = true;
            }
        }
        sum = 0;
    }
    bool mistake = true;
    while (mistake)
    {
        generator = std::uniform_real_distribution<double>(0.0, sumFit);
        double rip2 = generator(engine);
        mistake = false;
        foundCorrect = false;
        while (!foundCorrect)
        {
            double sum = 0;
            index2++;
            for (int j = 0; j < index2; j++)
            {
                sum += fitness[j];
                if (rip2 <= sum)
                {
                    foundCorrect = true;
                    if (index2 == index1)
                    {
                        mistake = true;
                        index2 = 0;
                    }
                }
            }
            sum = 0;
        }
    }
    parent = std::make_pair(index1 - 1, index2 - 1);
    return parent;
}

std::vector < std::pair<double, double>> RIP(std::vector<std::pair<double, double>> populationBeforeSelection)
{
    std::vector<double> fit;
    double sumFit = 0;
    for (int i = 0; i < 6; i++)
    {
        fit.push_back(2 - fitness(populationBeforeSelection[i]));
        sumFit += 2 - fitness(populationBeforeSelection[i]);
    }
    std::mt19937 engine(std::random_device{}());
    auto generator = std::uniform_real_distribution<double>(0.0, sumFit);
    double rip1 = generator(engine);
    int index1 = 0, index2 = 0;
    bool foundCorrect = false;
    while (!foundCorrect)
    {
        double sum = 0;
        index1++;
        for (int j = 0; j < index1; j++)
        {
            sum += fit[j];
            if (rip1 <= sum)
            {
                foundCorrect = true;
            }
        }
        sum = 0;
    }
    bool mistake = true;
    while (mistake)
    {
        generator = std::uniform_real_distribution<double>(0.0, sumFit);
        double rip2 = generator(engine);
        mistake = false;
        foundCorrect = false;
        while (!foundCorrect)
        {
            double sum = 0;
            index2++;
            for (int j = 0; j < index2; j++)
            {
                sum += fit[j];
                if (rip2 <= sum)
                {
                    foundCorrect = true;
                    if (index2 == index1)
                    {
                        mistake = true;
                        index2 = 0;
                    }
                }
            }
            sum = 0;
        }
    }
    std::vector<std::pair<double, double>> populationAfterSelection;
    for (int i = 0; i < 6; i++)
    {
        if ((i != index1 - 1) && (i != index2 - 1))
        {
            populationAfterSelection.push_back(populationBeforeSelection[i]);
        }
    }
    return populationAfterSelection;
}

std::vector<std::pair<double, double>> iteration(std::vector<std::pair<double, double>> unitVector)
{
    std::vector<std::pair<double, double>> newVector = unitVector;
    std::mt19937 engine(std::random_device{}());
    std::vector<double> fit(4);
    double maxFit = 0;
    double sumFit = 0;
    for (int i = 0; i < 4; i++)
    {
        fit[i] = fitness(unitVector[i]);
        sumFit += fitness(unitVector[i]);
    }
    std::pair<int, int> parent = parentGenerator(fit, sumFit);
    std::pair<std::pair<double, double>, std::pair<double, double>> newChildren = children(unitVector[parent.first], unitVector[parent.second]);
    newVector.push_back(newChildren.first);
    newVector.push_back(newChildren.second);
    for (int i = 0; i < 6; i++)
    {
        auto generator = std::uniform_real_distribution<double>(0.0, 1.0);
        double mutationChance = generator(engine);
        if (mutationChance < 0.25)
        {
            newVector[i] = mutation(newVector[i]);
        }
    }
    std::vector<std::pair<double, double>> vectorAfterDelete = RIP(newVector);
    sumFit = 0;
    unitVector = vectorAfterDelete;
    for (int i = 0; i < 4; i++)
    {
        fit[i] = fitness(unitVector[i]);
        sumFit += fitness(unitVector[i]);
    }
    for (int i = 0; i < 4; i++)
    {
        maxFit = std::max(maxFit, fit[i]);
    }
    for (int i = 0; i < 4; i++)
    {
        out << "\t\t|" << std::fixed << std::setprecision(5) << unitVector[i].first << " \t|" << unitVector[i].second << " \t|" << fit[i] << " \t|" << std::setw(16) << "|\n";
    }
    out << "\t\t|" << std::setw(16) << "|" << std::setw(16) << "|" << std::setw(16) << "|" <<  maxFit << std::setw(8) << "|" << sumFit / 4 << "\n";
    out << "---------------------------------------------------------------------------------------\n";
    return unitVector;
}

void geneticsAlgorithm(int generations)
{
    out << "Generation" << "\t| X\t\t| Y\t\t| FIT\t\t| Max" << std::setw(17) << "|Middl\n";
    out << "---------------------------------------------------------------------------------------\n";
    std::mt19937 engine(std::random_device{}());
    std::vector<std::pair<double, double>> unitVector;
    unitVector.resize(4);
    unitVector = Generator();
    for (int i = 1; i <= generations; i++)
    {
        out << i;
        unitVector = iteration(unitVector);
    }
}