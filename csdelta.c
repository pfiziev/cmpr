#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "utils.c"


void get_prob_vector(int bin_idx,
                     int * segmentations,
                     int n_datasets,
                     ldouble * metric,
                     int n_states,
                     int chrom_length,
                     ldouble * prob_vectors) {

    ldouble total = 0;

    for (int state_idx = 0; state_idx < n_states; state_idx ++) {

        // reset the probs vector
        prob_vectors[state_idx] = 0;

        // add up the probs for each dataset from its metric
        for (int dataset_idx = 0; dataset_idx < n_datasets; dataset_idx ++) {
            int dataset_state =  *(segmentations + dataset_idx * chrom_length + bin_idx);
            prob_vectors[state_idx] += *(metric + dataset_idx * n_states * n_states + dataset_state * n_states + state_idx);

//            printf("state= %d\td= %d\tds= %d\tp(y|x)= %lf\n", state_idx, dataset_idx, dataset_state,
//            *(metric + dataset_idx * n_states * n_states + dataset_state * n_states + state_idx));
        }

        total += prob_vectors[state_idx];
    }

    // normalize prob_vectors to sum to 1
    for (int state_idx = 0; state_idx < n_states; state_idx ++) {
        prob_vectors[state_idx] /= total;
    }
}


void get_overall_and_per_state_diff_score_for_multiple_groups(int * group_segmentations,
                                                              int n_groups,
                                                              int * group_offset,
                                                              int * datasets_per_group,

                                                              int chrom_length,

                                                              ldouble * metric,
                                                              int n_states,

                                                              ldouble * consistent_state_probs,

                                                              int compute_per_state_scores,
                                                              int compute_differential_from_closest_group,

                                                              // results variables
                                                              ldouble * overall_scores,
                                                              ldouble * per_state_scores,
                                                              ldouble * constitutive_scores,
                                                              int * closest_state_per_group,
                                                              int * closest_constitutive_state
                                                              ) {
    time_t timer;
    char buffer[26];
    struct tm* tm_info;

//    print_array_int(group_offset, n_groups);
//    print_array_int(datasets_per_group, n_groups);

    ldouble ** group_to_group_distance = matrix(n_groups, n_groups);

    ldouble * all_summed_probs = new_array(n_states);
    ldouble * rest_summed_probs = new_array(n_states);

    ldouble ** prob_vectors = matrix(n_groups, n_states);

    for (int bin_idx = 0; bin_idx < chrom_length; bin_idx ++) {
//    for (int bin_idx = 0; bin_idx < 1; bin_idx ++) {

        if (! (bin_idx % 20000)) {

            time(&timer);
            tm_info = localtime(&timer);

            strftime(buffer, 26, "%Y-%m-%d %H:%M:%S", tm_info);

            printf("[%s] Processed bins: %d out of %d\n", buffer, bin_idx, chrom_length);
        }

        for (int state_idx = 0; state_idx < n_states; state_idx ++ ) {
            all_summed_probs[state_idx] = 0;
        }


        for (int g_idx = 0; g_idx < n_groups; g_idx ++) {

            get_prob_vector(bin_idx,
                            group_segmentations + group_offset[g_idx] * chrom_length,
                            datasets_per_group[g_idx],
                            metric + group_offset[g_idx] * n_states * n_states,
                            n_states,
                            chrom_length,
                            prob_vectors[g_idx]);

//            printf("Group: %d\n", g_idx);
//            print_array(prob_vectors[g_idx], n_states);

            for (int state_idx = 0; state_idx < n_states; state_idx ++ ) {
                all_summed_probs[state_idx] += prob_vectors[g_idx][state_idx];

                // store the sum of all probs as constitutive score for this state
                *(constitutive_scores + state_idx * chrom_length + bin_idx) += prob_vectors[g_idx][state_idx];
            }
        }

        // find the state with the highest constitutive score and store it as the best guess for closest_constitutive_state
        ldouble best_constitutive_score = 0;

        for (int state_idx = 0; state_idx < n_states; state_idx ++ ) {

            ldouble constitutive_score = *(constitutive_scores + state_idx * chrom_length + bin_idx);

            if (constitutive_score > best_constitutive_score) {
                best_constitutive_score = constitutive_score;
                closest_constitutive_state[bin_idx] = state_idx;
            }
        }

        // compute the group to group distance
        if (compute_differential_from_closest_group) {
            for (int g_idx_1 = 0; g_idx_1 < n_groups; g_idx_1 ++) {
                for (int g_idx_2 = g_idx_1 + 1; g_idx_2 < n_groups; g_idx_2 ++) {

                    ldouble sym_kl_div = symmetric_KL_divergence(prob_vectors[g_idx_1],
                                                                 prob_vectors[g_idx_2],
                                                                 n_states);

                    group_to_group_distance[g_idx_1][g_idx_2] = sym_kl_div;
                    group_to_group_distance[g_idx_2][g_idx_1] = sym_kl_div;
                }
            }
        }

        for (int g_idx_1 = 0; g_idx_1 < n_groups; g_idx_1 ++ ) {
            ldouble * g_prob = prob_vectors[g_idx_1];

            if (compute_differential_from_closest_group) {
                // find the closest group
                int closest_group_idx = -1;
                ldouble min_dist = -1;

                for (int g_idx_2 = 0; g_idx_2 < n_groups; g_idx_2 ++ ) {
                    if (g_idx_2 == g_idx_1) {
                        continue;
                    }

                    if (min_dist == -1 || group_to_group_distance[g_idx_1][g_idx_2] < min_dist) {
                        min_dist = group_to_group_distance[g_idx_1][g_idx_2];
                        closest_group_idx = g_idx_2;
                    }
                }

                *(overall_scores + g_idx_1 * chrom_length + bin_idx) = min_dist;

                for (int state_idx = 0; state_idx < n_states; state_idx ++ ) {
                    *(per_state_scores + g_idx_1 * n_states * chrom_length + state_idx * chrom_length + bin_idx) = g_prob[state_idx] - prob_vectors[closest_group_idx][state_idx];
                }
            } else {

                for (int state_idx = 0; state_idx < n_states; state_idx ++ ) {
                    rest_summed_probs[state_idx] = (all_summed_probs[state_idx] - g_prob[state_idx]) / (n_groups - 1);
                }

                *(overall_scores + g_idx_1 * chrom_length + bin_idx) = symmetric_KL_divergence(g_prob, rest_summed_probs, n_states);

                for (int state_idx = 0; state_idx < n_states; state_idx ++ ) {
                    *(per_state_scores + g_idx_1 * n_states * chrom_length + state_idx * chrom_length + bin_idx) = g_prob[state_idx] - rest_summed_probs[state_idx];
                }

            }

            // find the best state for this group
            ldouble best_score = -1;
            int * cur_group_segmentations = group_segmentations + group_offset[g_idx_1] * chrom_length;

            for (int dataset_idx = 0; dataset_idx < datasets_per_group[g_idx_1]; dataset_idx ++ ) {

                int state_idx = *(cur_group_segmentations + dataset_idx * chrom_length + bin_idx);
                ldouble c_score = KL(g_prob, consistent_state_probs + g_idx_1 * n_states * n_states + state_idx * n_states, n_states);

                if (best_score == -1 || c_score < best_score) {
                    best_score = c_score;
                    *(closest_state_per_group + g_idx_1 * chrom_length + bin_idx) = state_idx;
                }
            }
        }



    }


    free_matrix(group_to_group_distance, n_groups);
    free_matrix(prob_vectors, n_groups);

    free(all_summed_probs);
    free(rest_summed_probs);
}
