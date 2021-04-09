#include "mpi.h"
#include <stdio.h>
#include "mmio.h"
#include "mmio.c"
#include <stdlib.h>
// global rank variable
void dispatchInfo(int max,int i, int j, int minRowsPerProcessor, int indexs[2]);
void calculateResults(int* rowInfo, int* colInfo, double* values, double* vectorValues, int num_rows, int buffer_size, int start_row, int start_col, double *reuslt);
int rank;
const int DIS_TARGET_DEF = 0;
const int DIS_TARGET_SYM = 1;
const int MASTER = 0;
int main(int argc, char* args[])

{
	//  mmio test
	MM_typecode matcode;
	int colNum, rowNum, nzerosNum = 0;
	MM_typecode matcodeVector;
	int colNumVector, rowNumVector, nzerosNumVector = 0;
	//int  pimatm = 0;
	int i, n;
	//int istage, irow, icol, jrow, jcol, iproc, jproc, index, Proc_Id ;
	//int A_Bloc_MatrixSize, B_Bloc_MatrixSize;
	//int NoofRows_A, NoofCols_A, NoofRows_B, NoofCols_B;
	//int NoofRows_BlocA, NoofCols_BlocA, NoofRows_BlocB, NoofCols_BlocB;
	//int Local_Index, Global_Row_Index, Global_Col_Index;
	//int Matrix_Size[4];
	//int source, destination, send_tag, recv_tag, Bcast_root;

	//float** Matrix_A, ** Matrix_B, ** Matrix_C;
	//float* A_Bloc_Matrix, * B_Bloc_Matrix, * C_Bloc_Matrix, * Temp_BufferA;

	//float* MatA_array, * MatB_array, * MatC_array;
	double** subMatrixes;
	double* verctor = NULL;
	int** rowIndexs;
	int** colIndexs;
	int* maxSizeArray;
	int* minIndexArray;
	int* maxIndexArray;
	int dispatchTarget[2];
	FILE* fp;
	FILE* fpv;

	int	 MatrixA_FileStatus = 1, VectorB_FileStatus = 1;

	int size;
	
	int row, col = 0;
	double value;

	MPI_Init(&argc, &args);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//symTargetProcessorsize = 4;

	subMatrixes = (double**)malloc(size * sizeof(double*));
	rowIndexs = (int**)malloc(size * sizeof(int*));
	colIndexs = (int**)malloc(size * sizeof(int*));
	maxSizeArray = (int*)calloc(size, sizeof(int));
	minIndexArray = (int*)calloc(size, sizeof(int));
	maxIndexArray = (int*)calloc(size, sizeof(int));

	// how many rows each processor can get
	int minRowsPerProcessor;
	int bufferSize;
	int maxProcessorId = size - 1;
	printf("init %d \r\n ", rank);

	// prepartion phase
	if (rank == MASTER)
	{
		printf("master read %d \r\n ", rank);
		if ((fp = fopen("sherman1.mtx", "r")) == NULL) {
			MatrixA_FileStatus = 0;
		}

		if (MatrixA_FileStatus != 0) {
			mm_read_banner(fp, &matcode);
			mm_read_mtx_crd_size(fp, &rowNum, &colNum, &nzerosNum);
			minRowsPerProcessor = rowNum / size;
			// the matrix is symatric, need to double the buffer.
			bufferSize = 4 * nzerosNum / size;

			// allocate memory
			for (int i = 0; i < size; i ++) {
				subMatrixes[i] = (double*)malloc(bufferSize * sizeof(double));
				rowIndexs[i] = (int*)malloc(bufferSize * sizeof(int));
				colIndexs[i] = (int*)malloc(bufferSize * sizeof(int));
			}

			for (i = 0; i < nzerosNum; i++)
			{
				fscanf(fp, "%d %d %lg\n", &row, &col, &value);
				row--;  /* adjust from 1-based to 0-based */
				col--;

				dispatchInfo(maxProcessorId ,row, col, minRowsPerProcessor, dispatchTarget);
				if (row == 424 || col == 424) {
					printf("=============%d %d %lg %d %d\n", row, col, value , dispatchTarget[0], dispatchTarget[1]);
				}
				int targetProcessor = dispatchTarget[DIS_TARGET_DEF];
				int current_index = maxSizeArray[targetProcessor];
				subMatrixes[targetProcessor][current_index] = value;
				rowIndexs[targetProcessor][current_index] = row;
				colIndexs[targetProcessor][current_index] = col;
				if (current_index == 0) {
					minIndexArray[targetProcessor] = col;
					maxIndexArray[targetProcessor] = col;
				}
				else {
					if (minIndexArray[targetProcessor] > col) {
						minIndexArray[targetProcessor] = col;
					}
					if (maxIndexArray[targetProcessor] < col) {
						maxIndexArray[targetProcessor] = col;
					}
				}

				current_index++;
				maxSizeArray[targetProcessor] = current_index;
				int symTargetProcessor = dispatchTarget[DIS_TARGET_SYM];
				if (symTargetProcessor!=-1) {
					current_index = maxSizeArray[symTargetProcessor];
					subMatrixes[symTargetProcessor][current_index] = value;
					rowIndexs[symTargetProcessor][current_index] = col;
					colIndexs[symTargetProcessor][current_index] = row;
					current_index++;
					maxSizeArray[symTargetProcessor] = current_index;
					if (minIndexArray[symTargetProcessor] > row) {
						minIndexArray[symTargetProcessor] = row;
					}
					if (maxIndexArray[symTargetProcessor] < row) {
						maxIndexArray[symTargetProcessor] = row;
					}
				}

			}
			fclose(fp);
			
			
		}
		if ((fpv = fopen("sherman1_rhs1.mtx", "r")) == NULL) {
			VectorB_FileStatus = 0;
		}
		if (VectorB_FileStatus) {
			mm_read_banner(fpv, &matcodeVector);
			mm_read_mtx_array_size(fpv, &nzerosNumVector, &rowNumVector);
			verctor = (double*)malloc(nzerosNumVector * sizeof(double));
			for (i = 0; i < nzerosNumVector; i++)
			{
				fscanf(fpv, "%lg\n", &verctor[i]);
			}
			for (i = 0; i < nzerosNumVector; i++) {
				//printf("%e\r\n", verctor[i]);
			}
		}
		else {
			exit(1);
		}
		printf("master read end %d \r\n ", rank);
		// debug
		//for (int i = 0; i < size; i++) {
			//printf("proceor %d , %d\r\n ", i, maxSizeArray[i]);
			//for (int j = 0;j < maxSizeArray[i]; j++) {
				//if (rowIndexs[i][j] >= colIndexs[i][j]) {
				//	printf("proceor %d,%d,%d\r\n", rowIndexs[i][j], colIndexs[i][j], subMatrixes[i][j]);
				//}
				//else {
				//	printf("proceor %d,%d,%d <-\r\n", rowIndexs[i][j], colIndexs[i][j], subMatrixes[i][j]);
				//}
			//}
		//}
	};
	printf("brocast %d \r\n ",rank);
	// brocast
	MPI_Bcast(&minRowsPerProcessor, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&rowNum, 1, MPI_INT, 0, MPI_COMM_WORLD);
	printf("brocast end  %d \r\n ", rank);
	MPI_Status status;
	int blockSize, minIndex, maxIndex;
	double *values;
	double* partical_vector_values;
	int *colInfo;
	int *rowInfo;
	if (rank == MASTER)
	{
		printf("MASTER send  %d \r\n ", size);
		for (int rankIndex = 1; rankIndex < size; rankIndex++) {
			int blockSize = maxSizeArray[rankIndex];
			int vectorSize = maxIndexArray[rankIndex] - minIndexArray[rankIndex] +1;
			printf("INITIAL!!: Processor %d send to %d the message with %d values. start at %d, end at %d ,length  %d \n", rank, rankIndex, blockSize, minIndexArray[rankIndex], maxIndexArray[rankIndex], vectorSize);
			// metadata:the number of entries sent to a worker processor
			MPI_Send(&maxSizeArray[rankIndex], 1, MPI_INT, rankIndex, NULL, MPI_COMM_WORLD);
			// metadata:the start of column index
			MPI_Send(&minIndexArray[rankIndex], 1, MPI_INT, rankIndex, NULL, MPI_COMM_WORLD);
			// metadata:the end of column index
			MPI_Send(&maxIndexArray[rankIndex], 1, MPI_INT, rankIndex, NULL, MPI_COMM_WORLD);
			// data: row info
			MPI_Send(rowIndexs[rankIndex], blockSize, MPI_INT, rankIndex, NULL, MPI_COMM_WORLD);
			//for (int j = 0; j < blockSize; j++) {
			//	printf(" %d \r\n ", rowIndexs[rankIndex][j]);
			//}
			// data: column info
			MPI_Send(colIndexs[rankIndex], blockSize, MPI_INT, rankIndex, NULL, MPI_COMM_WORLD);
			// data:values
			MPI_Send(subMatrixes[rankIndex], blockSize, MPI_DOUBLE, rankIndex, NULL, MPI_COMM_WORLD);
			int startIndex = minIndexArray[rankIndex];
			printf("===startIndex  %d \r\n ", startIndex);
			MPI_Ssend(&verctor[startIndex], vectorSize, MPI_DOUBLE, rankIndex, NULL, MPI_COMM_WORLD);
			
		}
		double* result = (double*)calloc(rowNum, sizeof(double));
		calculateResults(rowIndexs[0], colIndexs[0], subMatrixes[0], verctor, minRowsPerProcessor, maxSizeArray[0], 0, minIndexArray[0], result);
		//for (int j = 0; j < minRowsPerProcessor; j++) {
		//	printf(" %d  %f \r\n ", j,result[j]);
		//}
		// receive from works
		for (int rankIndex = 1; rankIndex < size; rankIndex++) {
			int offset = rankIndex * minRowsPerProcessor;
			int num_rows = minRowsPerProcessor;
			int start_row = rankIndex * minRowsPerProcessor;
			if (rankIndex == size - 1) {
				num_rows = rowNum - start_row;
			}
			MPI_Recv(&result[offset], num_rows, MPI_DOUBLE, rankIndex, NULL, MPI_COMM_WORLD, &status);
		}
		FILE* outf = NULL;

		outf = fopen("result.txt", "w+");
		fprintf(outf, "This is result\n");
		for (int j = 0; j < rowNum; j++) {
			fprintf(outf, "%d %e\n", j,result[j]);
		}

	}
	else {
		MPI_Recv(&blockSize, 1, MPI_INT, 0, NULL, MPI_COMM_WORLD, &status);
		MPI_Recv(&minIndex, 1, MPI_INT, 0, NULL, MPI_COMM_WORLD, &status);
		MPI_Recv(&maxIndex, 1, MPI_INT, 0, NULL, MPI_COMM_WORLD, &status);
		int vectorSize = maxIndex - minIndex+1;

		printf("INITIAL!!: Processor %d has received the message with %d values. min %d, max %d, vector length %d \n", rank, blockSize, minIndex, maxIndex, vectorSize);
		//printf("INITIAL!!: Processor %d has received the message with %d values.\n", rank, 1);
		//blockSize = 1;
		rowInfo = (int*)malloc(sizeof(int) * blockSize);
		colInfo = (int*)malloc(blockSize * sizeof(int));
		values = (double*)malloc(blockSize * sizeof(double));
		partical_vector_values = (double*)malloc(vectorSize * sizeof(double));
		MPI_Recv(rowInfo, blockSize, MPI_INT, 0, NULL, MPI_COMM_WORLD, &status);
		MPI_Recv(colInfo, blockSize, MPI_INT, 0, NULL, MPI_COMM_WORLD, &status);
		MPI_Recv(values, blockSize, MPI_DOUBLE, 0, NULL, MPI_COMM_WORLD, &status);
		MPI_Recv(partical_vector_values, vectorSize, MPI_DOUBLE, 0, NULL, MPI_COMM_WORLD, &status);
		for (int j = 0; j < blockSize; j ++) {
			//printf("->%d, %d, %d, %e \r\n ", rank ,rowInfo[j], colInfo[j], values[j]);
			//printf(" %d, %d \r\n ", rowInfo[j], colInfo[j]);
		//	printf(" %d \r\n ", rowInfo[j]);
		}
		for (int j = 0; j < vectorSize; j++) {
			if (rank == 1) {
				//printf("=V %d= %d  %f \r\n ", rank, j + minIndex, partical_vector_values[j]);
			}
		}
		int num_rows = minRowsPerProcessor;
		int start_row = rank * minRowsPerProcessor;
		if (rank == size - 1) {
			num_rows = rowNum - start_row;
		}
		double* result = (double*)calloc(num_rows, sizeof(double));
		calculateResults(rowInfo, colInfo, values, partical_vector_values, num_rows, blockSize, start_row, minIndex, result);
		for (int j = 0; j < num_rows; j++) {
			//printf("=%d= %d  %e  \r\n ", rank, j + start_row, result[j]);
			//printf("=== %d  %f \r\n ", j, partical_vector_values[j]);
		}
		MPI_Send(result, num_rows, MPI_DOUBLE, 0, NULL, MPI_COMM_WORLD);
	}



	MPI_Finalize();
	return  0;
}

void calculateResults(int *rowInfo, int*colInfo, double *values, double *vectorValues, int num_rows, int buffer_size,int start_row, int start_col,double *result) {
	printf("calculateResults!!: num_rows %d ,buffer_size %d ,start_row %d ,start_col %d ,.\n", num_rows, buffer_size, start_row, start_col);
	//
	for (int i = 0; i < buffer_size; i ++){
		
		int row = rowInfo[i];
		int col = colInfo[i];
		double value = values[i];
		
		int new_row = row - start_row;
		int new_col = col - start_col;
		result[new_row] = result[new_row] + value * vectorValues[new_col];
		
		if (row == 424) {
			printf("[%d] row %d ,col %d ,value %f,vector[%d] %e result[%d] %e \n", rank,row, col, value, new_col,vectorValues[new_col], new_row,result[new_row]);
		}
		if (col == 424) {
			//printf("**colInfo[i] %d ,start_row %d ,num_rows %d , row %d , col %d \n", colInfo[i], start_row, num_rows, row, col);
		}
		if (colInfo[i] >= start_row && colInfo[i] < start_row + num_rows && row!= col) {
			new_row = col - start_row;
			new_col = row - start_col;
			result[new_row] = result[new_row] + value * vectorValues[new_col];
			if (col == 424) {
				printf(">[%d]  row %d ,col %d ,value %e , vector[%d] %e , result[%d] %e \n", rank, col, row, value, new_col,vectorValues[new_col], new_row,result[new_row]);
			}
		}
	}
	return result;
}

void dispatchInfo(int max,int i, int j, int minRowsPerProcessor, int indexs[2]) {
	int proId = i / minRowsPerProcessor;
	proId = proId > max ? max : proId;
	indexs[DIS_TARGET_DEF] = proId;
	int symmetric_proId = j / minRowsPerProcessor;
	symmetric_proId = symmetric_proId > max ? max : symmetric_proId;
	if (indexs[DIS_TARGET_DEF] == symmetric_proId) {
		indexs[DIS_TARGET_SYM] = -1;
	}
	else {
		indexs[DIS_TARGET_SYM] = symmetric_proId;
	}
}