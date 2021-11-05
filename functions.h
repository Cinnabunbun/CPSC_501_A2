#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>

void scaleSamples(std::vector<double> *samples, int arraySize);

void convolve(double *convolutionArray, std::vector<double> *inputSamples, std::vector<double> *IRSamples, int inputSize, int IRSize);

std::vector<double>* readWavFile(int channels, char *filename);

void readWavFileHeader(int channels, char *inputFilename);

void writeWavFile(double *outputArray, int outputArraySize, int channels, char *filename);

void writeWavFileHeader(int channels, int numberSamples, double outputRate, FILE *outputFile);

size_t fwriteIntLSB(int data, FILE *stream);
int freadIntLSB(FILE *stream);
size_t fwriteShortLSB(short int data, FILE *stream);
short int freadShortLSB(FILE *stream);

#endif