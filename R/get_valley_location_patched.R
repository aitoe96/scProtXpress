###############################################################################
# get_valley_location_patched.R
#
# Local copy of ADTnorm::get_valley_location() (ADTnorm v1.0,
# https://github.com/yezhengSTAT/ADTnorm, R/get_valley_location.R), re-exposed
# under the internal name `.get_valley_location_patched()`.
# `ADTnorm::get_valley_location()` is NOT modified.
#
# (a) Why
#     When a sample has a single peak and neither a valley nor a shoulder is
#     found to its right, ADTnorm::get_valley_location() calls
#     stop("No valley is detected ...") and the whole marker is lost for
#     every sample. With cell-type balanced detection inputs this happens
#     more often (scProtXpress 0.3.0).
#
# (b) Changes relative to ADTnorm::get_valley_location() (v1.0)
#     - That stop() is replaced by a missing valley for the sample (NA),
#       which ADTnorm::landmark_fill_na() imputes, as it already does for
#       samples without peaks. The samples are returned in the new element
#       `no_valley_samples` of the output list.
#     - In the `multi_sample_per_batch` branch, the line
#         cell_x_feature %>% dplyr::filter(batch == batch_each) %$% sample %>% unique
#       is written in base R with the same result
#         unique(cell_x_feature$sample[which(cell_x_feature$batch == batch_each)])
#       because scProtXpress imports neither dplyr nor magrittr's %$%.
#     Nothing else is changed (roxygen header of the original omitted).
#
# License: GPL-3 (inherited from ADTnorm; see ADTnorm LICENSE.md).
# Original author: Ye Zheng <yzheng23@fredhutch.org>.
###############################################################################

#' @noRd
.get_valley_location_patched = function(cell_x_adt = NULL, cell_x_feature = NULL, adt_marker_select = NULL, peak_mode_res = NULL, shoulder_valley = TRUE, positive_peak = NULL, multi_sample_per_batch = FALSE, adjust = 1.5, min_fc = 20, shoulder_valley_slope = -1, lower_peak_thres = 0.01, neg_candidate_thres = asinh(10/5 + 1), arcsinh_transform_flag = TRUE) {

    peak_landmark_list = peak_mode_res
    no_valley_samples = character(0) ## scProtXpress: samples whose valley is left NA (see header)

    ## tag sample if it is the only sample within each batch that have different number of peaks.
    # batch_num = cell_x_feature$batch %>% unique %>% length
    tag_row = vector("list", length = nrow(peak_landmark_list))
    names(tag_row) = rownames(peak_landmark_list)

    if(multi_sample_per_batch && ncol(peak_landmark_list) > 1){ ## multiple sample per batch remove the positive peak that only appear in one sample
        for(batch_each in unique(cell_x_feature$batch)){
            sample_index = unique(cell_x_feature$sample[which(cell_x_feature$batch == batch_each)]) ## scProtXpress: base R (see header)
            ## if there is only one second peak, remove it
            if(length(sample_index) > 2){
                col_index = which(colSums(!is.na(peak_landmark_list[sample_index, ])) == 1) ## which col has only 1 peak per batch
                # common_peak_col = which(colSums(!is.na(peak_landmark_list[sample_index, ])) != 1) #ncol(peak_landmark_list[sample_index, ]) - length(col_index)
                for(c in col_index){ ## across those columns that only have 1 peak per batch
                    row_index = sample_index[which(!is.na(peak_landmark_list[sample_index, c]))] ## get the sample name
                    tag_row[[row_index]] = c(tag_row[[row_index]], c) #common_peak_n #non-zero value, record the peak number that is common across samples in this batch
                }
            }
        }
    }

    valley_location_list = matrix(NA, nrow = nrow(peak_landmark_list), ncol = max(1, ncol(peak_landmark_list) - 1))
    rownames(valley_location_list) = cell_x_feature$sample %>% levels()
    for (sample_name in cell_x_feature$sample %>% levels()) {
        # cell_ind = which(cell_x_feature$sample == sample_name)
        cell_ind_tmp = which(cell_x_feature$sample == sample_name)
        cell_notNA = which(!is.na(cell_x_adt[cell_ind_tmp, adt_marker_select]))
        cell_ind = cell_ind_tmp[cell_notNA]
        if(length(cell_ind) > 0){
            peak_landmark = peak_landmark_list[sample_name, ]
            if(arcsinh_transform_flag){
                zero_prop = sum(cell_x_adt[cell_ind, adt_marker_select] < 2) / length(cell_x_adt[cell_ind, adt_marker_select])
            }else{
                zero_prop = sum(cell_x_adt[cell_ind, adt_marker_select] < neg_candidate_thres) / length(cell_x_adt[cell_ind, adt_marker_select])
            }

            ## check if user define single peak to be positive peak
            pos_marker_index = which(paste0("tmpName", positive_peak$ADT_index) == adt_marker_select)
            pos_sample_index = which(positive_peak$sample == sample_name)
            if (length(intersect(pos_marker_index, pos_sample_index)) > 0) {
                lower_valley = TRUE
            } else {
                lower_valley = FALSE
            }


            density_res = stats::density(
                cell_x_adt[which(cell_x_feature$sample == sample_name), adt_marker_select],
                adjust = adjust, na.rm = TRUE
            )
            x = density_res$x
            y = density_res$y

            sign_diff = sign(diff(y))
            diff_sign_diff = diff(sign_diff)
            peak = which(diff_sign_diff == -2) + 1
            valley = which(diff_sign_diff == 2) + 1

            x_valley = x[valley]
            y_valley = y[valley]
            x_peak = x[peak]
            real_peak = peak_landmark[!is.na(peak_landmark)] # peak
            np = length(real_peak)

            ## if this sample is tagged, choose the highest peak
            if (!is.null(tag_row[[sample_name]])) {## there are peaks need to be removed
                peak_landmark_y = c()
                for (peak_landmark_each in peak_landmark) {
                    peak_landmark_y = c(peak_landmark_y, y[which.min(abs(x - peak_landmark_each))]) ## get the density value for the peak location
                }
                common_peak_n = length(peak_landmark_y) - length(tag_row[[sample_name]])
                real_peak = real_peak[sort(order(peak_landmark_y, decreasing = T)[1:common_peak_n])] ## get the highest peak
                np = length(real_peak) ## update number of peak i.e, np
                peak_landmark_list[sample_name, ] = NA
                peak_landmark_list[sample_name, (1:ncol(peak_landmark_list))[-tag_row[[sample_name]]]] = real_peak ## remove the outlier peak
            }


            if (np > 1) { ## two peaks or more
                real_valley = c()
                for (i in 1:(np - 1)) {
                    tmp_valley = x_valley[(x_valley > real_peak[i]) & (x_valley < real_peak[i + 1])]
                    tmp_real_valley = tmp_valley[which.min(y_valley[(x_valley > real_peak[i]) & (x_valley < real_peak[i + 1])])]
                    if(length(tmp_real_valley) == 0){ ## no minimal point between two peak, ensure tmp_real_valley is not empty
                        tmp_real_valley = (real_peak[i] + real_peak[i+1]) / 2
                    }
                    if(i == 1 && shoulder_valley){
                        shoulder_cand_index = which(diff(y)/diff(x) > shoulder_valley_slope)
                        first_peak_index = (which(x > max(x_peak[1], real_peak[1])) %>% min) + 50
                        x_shoulder = x[shoulder_cand_index[shoulder_cand_index > first_peak_index][1]]
                        real_valley = c(real_valley, min(x_shoulder, tmp_real_valley, na.rm = T))
                    }else{
                        if(length(tmp_valley) == 0){ ## there is no local minimum that fall within two peaks
                            real_valley = c(real_valley, (real_peak[i] + real_peak[i+1]) / 2) ## use the midpoint of two peak location --- maybe shoulder peak instead?
                        }else{
                            real_valley = c(real_valley, tmp_real_valley) #c(real_valley, tmp_valley[which.min(y_valley[(x_valley > real_peak[i]) & (x_valley < real_peak[i + 1])])]) ## local minimum based on the density
                        }

                    }
                }
            } else if(np == 1) { ## one peak
                if (lower_valley == FALSE) { ## one peak is negative peak
                    real_valley =  x_valley[x_valley > real_peak[1] + 0.1][1]  #x[which(y < max(y) / min_fc)[which(y < max(y) / min_fc) > max(which(y == max(y)), which(x > real_peak[1]) %>% min())] %>% min()]
                    if(shoulder_valley){
                        ## if one peak & go for shoulder threshold
                        shoulder_cand_index = which(diff(y)/diff(x) > shoulder_valley_slope)
                        first_peak_index = (which(x > max(x_peak[1], real_peak[1])) %>% min) + 50
                        x_shoulder = x[shoulder_cand_index[shoulder_cand_index > first_peak_index][1]]
                        if(is.na(x_shoulder) & is.na(real_valley)){
                            ## scProtXpress: valley left NA (imputed by landmark_fill_na) instead of
                            ## stop("No valley is detected. Please consider increasing 'valley_density_adjust'.")
                            real_valley = NA_real_
                            no_valley_samples = c(no_valley_samples, sample_name)
                        }else{
                            real_valley = min(x_shoulder, real_valley, na.rm = T)
                        }
                        

                    }else{
                        ## check if no valley is detected due to shoulder peak
                        # if(length(y_valley[x_valley > real_peak[1]]) == 0 || (y_valley[x_valley >  real_peak[1] + 0.1][1] < 0.05)){
                        #     ## if one peak consider the shoulder point
                        #     shoulder_cand_index = which(diff(y)/diff(x) > shoulder_valley_slope)
                        #     first_peak_index = (which(x > max(x_peak[1], real_peak[1])) %>% min) + 50
                        #     x_shoulder = x[shoulder_cand_index[shoulder_cand_index > first_peak_index][1]]
                        #     real_valley = min(x_shoulder, real_valley, na.rm = T)
                        # }
                        if (length(y_valley[!is.na(x_valley) & x_valley > real_peak[1]]) == 0 || (y_valley[!is.na(x_valley) & x_valley > real_peak[1] + 0.1][1] < 0.05)) {
                            shoulder_cand_index = which(diff(y)/diff(x) > shoulder_valley_slope)
                            first_peak_index = (which(x > max(x_peak[1], real_peak[1])) %>% min) + 50
                            x_shoulder = x[shoulder_cand_index[shoulder_cand_index > first_peak_index][1]]
                            real_valley = min(x_shoulder, real_valley, na.rm = TRUE)
                        }
                    }
                    if (zero_prop > 0.8) {
                        real_valley = max(neg_candidate_thres, real_valley)
                    }
                } else { ## one peak is positive peak
                    real_valley = x[which((y < max(y) / min_fc) | (y < lower_peak_thres))[which((y < max(y) / min_fc) | (y < lower_peak_thres))< min(which(y == max(y)), which(x < real_peak[1]) %>% max())] %>% max()]
                    # peak_landmark_list[sample_name, ] = asinh(0/5 + 1)
                    # peak_landmark_list[sample_name, ncol(peak_landmark_list)] = real_peak[1]
                    # real_valley = min(real_valley, min(cell_x_adt[cell_ind, adt_marker_select])) #min(real_valley, asinh(0/5 + 1))
                }
            } else { ## no peak
                real_valley = NA
            }
            ## check if no valley is detected due to shoulder peak
            if (length(real_valley) == 0 || all(is.na(real_valley))) {
                if(length(real_peak) >= 2){ ## midpoint of two peak
                    real_valley = (real_peak[1] + real_peak[2]) / 2
                }else{## shoulder valley
                    shoulder_cand_index = which(diff(y)/diff(x) > shoulder_valley_slope)
                    first_peak_index = (which(x > max(x_peak[1], real_peak[1])) %>% min) + 50
                    x_shoulder = x[shoulder_cand_index[shoulder_cand_index > first_peak_index][1]]
                    if(!is.na(x_shoulder)){
                        real_valley = min(x_shoulder, real_valley, na.rm = T)
                    }
                }

            }
            valley_location_list[sample_name, 1:length(real_valley)] = real_valley
        }


    }

    ## update peak_landmark_list if needed
    if(any(colSums(!is.na(peak_landmark_list)) == 0)){ ## remove column that only has NA values
        col_index = which(colSums(!is.na(peak_landmark_list)) == 0)
        peak_landmark_update = matrix(NA, nrow(peak_landmark_list), ncol(peak_landmark_list) - length(col_index))
        i = 1
        for (c in 1:ncol(peak_landmark_list)) {
            if (!(c %in% col_index)) {
            peak_landmark_update[, i] = peak_landmark_list[, c]
            i = i + 1
            }
        }
        rownames(peak_landmark_update) = rownames(peak_landmark_list)
        peak_landmark_list = peak_landmark_update
    }
    return(list(valley_location_list = valley_location_list, peak_landmark_list = peak_landmark_list,
                no_valley_samples = no_valley_samples))
}
