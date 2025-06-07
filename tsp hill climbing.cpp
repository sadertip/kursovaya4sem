// tsp hill climbing.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <ctime>
#include <algorithm>
#include <random>
#include <numeric>  // ��������� ��� std::iota

void RandomfillMatrix(std::vector<std::vector<int>>& matrix) {
    int N = matrix.size();
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            matrix[i][j] = rand() % (50 - 1 + 1) + 2;  // ������� �����������
            matrix[j][i] = matrix[i][j];                    // ������ ����������� (���������)
        }
        matrix[i][i] = 0; 
    }
}

void FillMatrix(std::vector<std::vector<int>>& matrix) {
    int N = matrix.size();
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = i + 1; j < N; ++j)
        {
            int r;
            std::cin >> r;
            if (r < 0)
            {
                RandomfillMatrix(matrix);
                return;

            }
            matrix[i][j] = r;
            matrix[j][i] = matrix[i][j];  // ���������
        }
    }
}
void PrintMatrix(std::vector<std::vector<int>>& matrix) {
    for (const auto& row : matrix) {
        for (int val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "___________________________";
}
void PrintVec(std::vector<int> vec){
    for (const auto& item : vec) {
        std::cout << item << " ";
    }
}
int LenCount(std::vector<int> route, std::vector<std::vector<int>>& matrix){
    int first_city = route[0];
    int start_city = first_city;
    int total = 0;
    for (auto it = route.begin() + 1; it != route.end(); ++it) {
        total+= matrix[start_city][*it];
        start_city = *it;
    }
    return total+ matrix[start_city][first_city];
}
void twoOptSwap(std::vector<int>& path, int i, int k) {
    while (i < k) {
        std::swap(path[i], path[k]);
        i++;
        k--;
    }
}
std::vector<int> createVector(int n) {
    std::vector<int> vec(n);

    // ��������� ������� �� 0 �� n-1
    std::iota(vec.begin(), vec.end(), 0);

    // ������������
    std::shuffle(vec.begin(), vec.end(), std::mt19937{ std::random_device{}() });

    return vec;
}

std::vector<int> hillClimbingTSP(int max_iterations, std::vector<std::vector<int>>& matrix) {
    int n = matrix.size();
    std::vector<int> currentPath = createVector(n);
    int currentDistance = LenCount(currentPath,matrix);

    for (int iter = 0; iter < max_iterations; iter++) {
        bool improved = false;
        std::vector<int> bestNeighbor = currentPath;
        int bestNeighborDistance = currentDistance;
        // ���������� ���� ������� (Swap ���� �������)
        for (int i = 1; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                std::vector<int> neighbor = currentPath;
                std::swap(neighbor[i], neighbor[j]); // ������ ��� ������
                int neighborDistance = LenCount(neighbor,matrix);

                // ���� ����� ������ �������, ��������� �� ����
                if (neighborDistance < bestNeighborDistance) {
                    bestNeighbor = neighbor;
                    bestNeighborDistance = neighborDistance;
                    improved = true;
                }
            }
        }
        // ��������� ���� ������� ����� 2-opt
        for (int i = 1; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                std::vector<int> neighbor = currentPath;
                twoOptSwap(neighbor, i, j);
                int neighborDistance = LenCount(neighbor, matrix);

                if (neighborDistance < bestNeighborDistance) {
                    bestNeighbor = neighbor;
                    bestNeighborDistance = neighborDistance;
                    improved = true;
                }
            }
        }

        // ���� ����� ���������, ��������� � ������� ������
        if (improved) {
            currentPath = bestNeighbor;
            currentDistance = bestNeighborDistance;
            PrintVec(bestNeighbor);
            std::cout << std::endl << "score " << bestNeighborDistance << std::endl;
        }
        else {
            break; // ��������� ������� ���������
        }
    }


    return currentPath;
}



int main()
{
    size_t n;
    std::cin >> n;
    std::vector<std::vector<int>> roads(n, std::vector<int>(n));
    
    FillMatrix(roads); //��� �������� �� ����������� ������������, �����?
    PrintMatrix(roads);

    std::vector<int> route = hillClimbingTSP(100, roads);
    PrintVec(route);
    int score = LenCount(route, roads);
    std::cout<< std::endl<< score;

}
