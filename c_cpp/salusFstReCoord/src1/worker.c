#include <math.h>  // floor

#include "common.h"

static inline void transCorrd(double *new_pos, double pos_x, double pos_y, int_least16_t fovC, int_least16_t fovR) {
	/*
	    : pos_x - The original x coordinates
	    : pos_y - The original y coordinates
	    : fovR - The row number of fov
	    : fovC - The col number of fov
	*/
	double new_x = (fovC - (CenterFOV_COL - 1)) * FOV_USED_WIDTH - pos_y - 1;
	double new_y = pos_x + (fovR - (CenterFOV_ROW - 1) - 1) * FOV_USED_HEIGHT;
	new_pos[0] = floor(new_x * 100) / 100;
	new_pos[1] = floor(new_y * 100) / 100;
}

static inline char *strncpy_no_colon(char *dest, const char *src, size_t n) {
	size_t j = 0;  // destination index
	for (size_t i = 0; i < n && src[i] != '\0'; i++) {
		if (src[i] != ':') dest[j++] = src[i];
	}
	// for ( ; j < n; j++) dest[j] = '\0';
	if (j < n) dest[j] = '\0';
	return dest;
}

void worker(int_least16_t worker_id) {
	workerArray_t *worker = &Parameters.workerArray[worker_id];
	uint64_t index = 0;
	regmatch_t matches[2];
	char **splitSets = worker->tokens;
	char buffer[MAXFQIDLEN + 1] = {0};  // 81 => [0,80]
	char fovRC[9] = {0};                // R123C567
	for (uint64_t index = 0; index < JOBITEMSIZE; index++) {
		fstBCdata_t *fstBCdata_p = &worker->input_array[index];
		if (fstBCdata_p->name[0] == 0) {
			break;
		}
		fstBCoutput_t *fstBCoutput_p = &worker->output_array[index];
		memcpy(buffer, (char *)fstBCdata_p->name, sizeof(fstBCdata_p->name));
		assert(buffer[sizeof(fstBCdata_p->name)] == '\0');
		printf("%d\t[%s] <--\n", index, buffer);
		int_least16_t RowCol[2] = {0};
		double oldXY[2] = {0};
		double newXY[2];
		char *delim = ":";
		if (regexec(&Parameters.regex, buffer, 2, matches, 0) != 0) {
			continue;
		}
		//__builtin_dump_struct(&matches[0], &printf);
		//__builtin_dump_struct(&matches[1], &printf);
		size_t relen = matches[0].rm_eo - matches[0].rm_so;
		strncpy_no_colon(fovRC, buffer + matches[0].rm_so, relen);
		// fovRC[relen-1] = '\0';
		assert(fovRC[sizeof(fovRC) - 1] == '\0');
		RowCol[0] = atoi(fovRC + 1);
		RowCol[1] = atoi(fovRC + 5);
		if (unlikely(matches[1].rm_so == -1)) {
			delim = "_";
		} else {
			delim = ":";
		}
		assert(relen == (delim[0] == ':' ? sizeof(fovRC) : sizeof(fovRC) - 1));
		printf("%d\t[%s], delim:[%s], fov[%s]\n", index, buffer, delim, fovRC);
		const char *theDelim = delim;
		char *saveptr = NULL;
		char *token = strtok_r(buffer, theDelim, &saveptr);
		printf("-f- [%zu] [%zu] [%s]\n", buffer, token, token);
		if (unlikely(token == NULL)) {
			printf("-b->\t[%s], delim:[%s]\n", buffer, theDelim);
			break;
		} else {
			int_least16_t idx = 0;
			for (idx = 0; likely((token = strtok_r(NULL, theDelim, &saveptr)) != NULL); idx++) {
				splitSets[idx] = token;
				printf("-s- %d:[%zu] [%s]\n", idx, splitSets[idx], token);
			}
			oldXY[1] = atof(splitSets[idx - 1]);
			oldXY[0] = atof(splitSets[idx - 2]);
			printf("[%s],[%s]\n", splitSets[idx - 2], splitSets[idx - 1]);
			printf("[%s]->[%s]=(%d,%d), X:%f Y:%f\n", buffer, fovRC, RowCol[0], RowCol[1], oldXY[0], oldXY[1]);
			for (idx = 0; idx < MAXDELIMITEMS; idx++) {
				printf("-t- %d:[%zu] [%s]\n", idx, splitSets[idx], splitSets[idx]);
			}
		}
	}
}
