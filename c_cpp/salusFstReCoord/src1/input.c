#include <math.h>   // floor
#include <stdio.h>  // fprintf
#include <zlib-ng.h>

#include "kseq.h"
KSEQ_INIT(gzFile, zng_gzread)
#include "common.h"

void fqReader_init(void) {
	int rc;
	const char *pattern = "R[0-9]{3}(:)?C[0-9]{3}";
	if ((rc = regcomp(&Parameters.regex, pattern, REG_EXTENDED)) != 0) {
		regerror(rc, &Parameters.regex, Parameters.buffer, PARAMETERS_BUFFER_SIZE);
		printf("[x]=%d=regcomp(\"%s\") failed with '%s'.\n", rc, pattern, Parameters.buffer);
		exit(EXIT_FAILURE);
	}
	Parameters.ksfp = zng_gzopen(Parameters.inFastqFilename, "r");
	if (unlikely(Parameters.ksfp == NULL)) {
		fprintf(stderr, "[x]gzopen error on opening [%s]: %s.\n", Parameters.inFastqFilename, strerror(errno));
		exit(1);
	}
	rc = zng_gzbuffer(Parameters.ksfp, GZBUFSIZE);
	if (unlikely(rc != 0)) {
		fprintf(stderr, "[x]gzbuffer error: %s.\n", strerror(errno));
		exit(1);
	}
	Parameters.kseq = kseq_init(Parameters.ksfp);
}

void fqReader_destroy(void) {
	regfree(&Parameters.regex);
	kseq_destroy(Parameters.kseq);
	zng_gzclose(Parameters.ksfp);
}

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

void fill_worker(int_least16_t worker_id) {
	// ###### Below is from GPT, not reviewed yet. ######
	kseq_t *seq = Parameters.kseq;
	workerArray_t *worker = &Parameters.workerArray[worker_id];
	uint64_t index = 0;
	regmatch_t matches[2];

	while (index < JOBITEMSIZE && kseq_read(seq) >= 0) {
		fstBCdata_t *fstBCdata_p = &worker->input_array[index];
		fstBCoutput_t *fstBCoutput_p = &worker->output_array[index];
		strncpy(&fstBCdata_p->name, seq->name.s, sizeof(fstBCdata_p->name));
		/* seq->comment.s is discarded */
		strncpy((char *)fstBCdata_p->seq, seq->seq.s, sizeof(fstBCdata_p->seq));
		strncpy((char *)fstBCdata_p->qual, seq->qual.s, sizeof(fstBCdata_p->qual));

		// Extract the necessary fields from the FASTQ read ID
		char *read_id = seq->name.s;
		char fov[9];  // R123C456

		// Parse FOV and position
		int_least16_t row = 0, col = 0;
		double pos_x = 0.0, pos_y = 0.0;
		double new_pos[2];
		char delim = ':';
		fprintf(stderr, "%s\n", read_id);

		if (regexec(&Parameters.regex, read_id, 2, matches, 0) == 0) {
			size_t relen = matches[0].rm_eo - matches[0].rm_so;
			strncpy_no_colon(fov, read_id + matches[0].rm_so, relen);
			fov[relen] = '\0';
			row = atoi(fov + 1);
			col = atoi(fov + 5);
			if (unlikely(matches[0].rm_so == -1)) {
				delim = '_';
			} else {
				delim = ':';
			}
			printf("%.*s\n", 1, &delim);
			char *saveptr;
			char *token = strtok_r(read_id, &delim, &saveptr);
			if (unlikely(token == NULL)) {
				break;
			}
			char **split_set = fstBCdata_p->tokens;
			int_least16_t idx = 0;
			while (unlikely(token == NULL)) {
				split_set[idx++] = token;
				token = strtok_r(NULL, &delim, &saveptr);
			}
			pos_y = atof(split_set[idx - 1]);
			pos_x = atof(split_set[idx - 2]);
			printf("[%s]->[%s]=(%d,%d), X:%f Y:%f\n", read_id, fov, row, col, pos_x, pos_y);
			if ((FOV_X_MIN < pos_x && pos_x <= FOV_X_MAX) && (FOV_Y_MIN < pos_y && pos_y <= FOV_Y_MAX)) {
				pos_x -= FOV_X_MIN;
				pos_y -= FOV_Y_MIN;
				// transCorrd(new_pos, pos_x, pos_y, row, col);
				printf("-->gX:%.2f gY:%.2f\n", new_pos[0], new_pos[1]);
				snprintf(fstBCoutput_p->SpatiaStr, sizeof(worker->output_array[index].SpatiaStr), "%s %lf %lf", read_id, new_pos[0], new_pos[1]);
			}
		}
		index++;
		worker->fqSliceID++;  // Increment the slice ID for the worker
	}
	// Update the worker's flag to indicate the processing is complete
	atomic_store(&worker->flag, 1);
}
