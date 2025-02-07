
function [label_table] = compute_correctness_group(label_table)
    label_table.correctness_group = strings(height(label_table), 1);
    
    % Define the new trial correctness classes based on the updated criteria
    for i = 1:height(label_table)
        if sum(label_table.encoding_correctness(i, :), 2) == 3
            label_table.correctness_group(i) = "Correct_Trial";
        elseif sum(label_table.encoding_correctness(i, :), 2) < 2
            label_table.correctness_group(i) = "Incorrect_Trial";
        else
            label_table.correctness_group(i) = "Other";
        end
    end
end