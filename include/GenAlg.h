#pragma once
#include <iostream>
#include <random>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <fstream>

double fitness(std::pair<double, double> unit); // Фитнес функция

std::vector<std::pair<double, double>> Generator(); // Генерация начальной популяции

std::pair<std::pair<double, double>,
std::pair<double, double>> children(std::pair<double, double> parent1, std::pair<double, double> parent2); // Кроссовер

std::pair<double, double> mutation(std::pair<double, double> unit); // Мутация

std::pair<int, int> parentGenerator(std::vector<double> fitness, double sumFit); // Выбор родителей

std::vector < std::pair<double, double>> RIP(std::vector<std::pair<double, double>> populationBeforeSelection); // Селекция

std::vector<std::pair<double, double>> iteration(std::vector<std::pair<double, double>> unitVector); // Итерация

void geneticsAlgorithm(int generations); // Генетический алгоритм
