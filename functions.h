#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <complex>


typedef std::complex<double> cd;

void scaleDownIFFTResult(std::vector<double> &vec);

std::vector<double> complexArrayToDoubleArray(std::vector<cd> complex);

std::vector<cd> multComplexArrays(std::vector<cd> &arr1, std::vector<cd> &arr2);

int findFFTInputExponent(int outputSize);

void padZeros(std::vector<double> &v, int newSize);

void fft(std::vector<cd> &a, std::vector<cd> &A, int log2n);

void ifft(std::vector<cd> &a, std::vector<cd> &A, int log2n);

void scaleSamples(std::vector<double> *samples, int arraySize);

std::vector<double>* readWavFile(int channels, char *filename);

void readWavFileHeader(int channels, char *inputFilename);

void writeWavFile(double *outputArray, int outputArraySize, int channels, char *filename);

void writeWavFileHeader(int channels, int numberSamples, double outputRate, FILE *outputFile);

size_t fwriteIntLSB(int data, FILE *stream);
int freadIntLSB(FILE *stream);
size_t fwriteShortLSB(short int data, FILE *stream);
short int freadShortLSB(FILE *stream);

#endif