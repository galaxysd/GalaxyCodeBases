#include <math.h>  // floor

#include "common.h"

FORCE_INLINE void transCorrd(double *restrict ChipXY, const double *restrict FovXY, const int_least16_t *restrict FovRowCol) {
	double new_x = (FovRowCol[1] - (CenterFOV_COL - 1)) * FOV_USED_WIDTH - FovXY[1] - 1;
	double new_y = FovXY[0] + (FovRowCol[0] - (CenterFOV_ROW - 1) - 1) * FOV_USED_HEIGHT;
	double epsilon = 0.75 * OUTPUTPRECISION; /* eps in (0.5,1)/100, thus 0.75/100 */
	ChipXY[0] = epsilon + (floor(new_x * RECIPROCALPRECISION) * OUTPUTPRECISION);
	ChipXY[1] = epsilon + (floor(new_y * RECIPROCALPRECISION) * OUTPUTPRECISION);
}

FORCE_INLINE char *strncpy_no_colon(char *restrict dest, const char *restrict src, size_t n) {
	size_t j = 0;  // destination index
	for (size_t i = 0; i < n && src[i] != '\0'; i++) {
		if (src[i] != ':') dest[j++] = src[i];
	}
	// for ( ; j < n; j++) dest[j] = '\0';
	if (j < n) dest[j] = '\0';
	return dest;
}

void worker(int_least16_t worker_id) {
	workerArray_t *worker = &Parameters.worksQuene[worker_id];
	regmatch_t matches[2];
	char **splitSets = worker->tokens;
	// char* readName = malloc(91);          // for testing CHARsCPYSTR
	char readName[MAXFQIDLEN + 1] = {0};  // 81 => [0,80]
	char fovRC[9] = {0};                  // R123C567
	double GIVENunZoomRatio = (double)Parameters.unZoomRatio;
	double unZoomRatio = 0;
	for (uint64_t index = 0; index < JOBITEMSIZE; index++) {
		fstBCdata_t *fstBCdata_p = &worker->jobDatArray[index];
		if (fstBCdata_p->name[0] == 0) {
			break;
		}
		// fstBCoutput_t *fstBCoutput_p = &worker->output_array[index];
		ARRAYcpySTR(readName, fstBCdata_p->name);
		assert(readName[sizeof(fstBCdata_p->name)] == '\0');
#ifndef RELEASE
		fprintf(stderr, "###### %llu\t[%s] %d <--\n", index, readName, readName[sizeof(fstBCdata_p->name)]);
#endif
		int_least16_t RowCol[2] = {0};
		double oldXY[2] = {0.0};
		double newXY[2] = {0.0};
		char *delim = ":";
		if (regexec(&Parameters.regex, readName, 2, matches, 0) != 0) {
			fstBCdata_p->name[0] = 0;
			continue;
		}
		//__builtin_dump_struct(&matches[0], &printf);
		//__builtin_dump_struct(&matches[1], &printf);
		size_t relen = matches[0].rm_eo - matches[0].rm_so;
		strncpy_no_colon(fovRC, readName + matches[0].rm_so, relen);
		// fovRC[relen-1] = '\0';
		assert(fovRC[sizeof(fovRC) - 1] == '\0');
		// memcpy(fstBCdata_p->RowCol, fovRC, relen-1);
		fstBCdata_p->fov_row = (uint8_t)(RowCol[0] = atoi(fovRC + 1));
		fstBCdata_p->fov_column = (uint8_t)(RowCol[1] = atoi(fovRC + 5));
		if (unlikely(matches[1].rm_so == -1)) {
			delim = "_";
			unZoomRatio = 1.0;
		} else {
			delim = ":";
			unZoomRatio = GIVENunZoomRatio;
		}
		assert(relen == (delim[0] == ':' ? sizeof(fovRC) : sizeof(fovRC) - 1));
		// printf("%llu\t[%s], delim:[%s], fov[%s]\n", index, readName, delim, fovRC);
		const char *theDelim = delim;
		char *saveptr = NULL;
		char *token = strtok_r(readName, theDelim, &saveptr);
		// printf("-f- [%zu] [%zu] [%s]\n", readName, token, token);
		if (unlikely(token == NULL)) {
			// printf("-b->\t[%s], delim:[%s]\n", readName, theDelim);
			fstBCdata_p->name[0] = 0;
			break;
		} else {
			int_least16_t idx = 0;
			for (idx = 0; likely(token != NULL); token = strtok_r(NULL, theDelim, &saveptr), idx++) {
				splitSets[idx] = token;
				// printf("-s- %d:[%zu] [%s]\n", idx, splitSets[idx], token);
			}
			oldXY[1] = atof(splitSets[idx - 1]);
			oldXY[0] = atof(splitSets[idx - 2]);
#ifdef DEBUG
			fprintf(stderr, "oldXY: [%s],[%s] as [%.2f %.2f]\n", splitSets[idx - 2], splitSets[idx - 1], oldXY[0], oldXY[1]);
			/* for (idx = 0; idx < MAXDELIMITEMS; idx++) {
			    printf("-t- %d:[%zu] [%s]\n", idx, (void *)splitSets[idx], splitSets[idx]);
			} */
#endif
		}
		if (unlikely(unZoomRatio != 1.0f)) {
			oldXY[0] /= unZoomRatio;
			oldXY[1] /= unZoomRatio;
		}
#ifdef DEBUG
		fprintf(stderr, "[%s]->[%s]=RC(%d,%d), X:%f Y:%f\n", readName, fovRC, RowCol[0], RowCol[1], oldXY[0], oldXY[1]);
#endif
		if ((FOV_X_MIN < oldXY[0] && oldXY[0] <= FOV_X_MAX) && (FOV_Y_MIN < oldXY[1] && oldXY[1] <= FOV_Y_MAX)) {
			oldXY[0] -= FOV_X_MIN;
			oldXY[1] -= FOV_Y_MIN;
			transCorrd(newXY, oldXY, RowCol);
			// printf("-->gX:%.2f gY:%.2f\n", newXY[0], newXY[1]);
			memcpy(fstBCdata_p->newXY, newXY, sizeof(newXY));
#ifdef DEBUG
			char readSeq[BARCODELEN + 1];
			ARRAYcpySTR(readSeq, fstBCdata_p->seq);
			fprintf(stderr, "->SpatiaStr:[%s %.2f %.2f]\n", readSeq, fstBCdata_p->newXY[0], fstBCdata_p->newXY[1]);
#endif
		} else {
			fstBCdata_p->name[0] = 0;
		}
	}
}
