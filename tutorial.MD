# algorithm.py

This file stores the main algorithm code.

The algorithm beings at the top of algorithm.py

We take the fasta file and create a list where each element is a string form of a row in the data

Next, we preprocess the data (preprocess function) which transforms anything that is in the problems list to a 'B' character

Then, we get only the unstable sites (sites where H(X)>=epsilon)
From this point forward, you will notice that we have an indices and cols list:
  This is because we want to be able to match columns to which site they represent in the original data
  So, the indices list looks like [472,576,...] where each element is the index of the column of the same element in the corresponding cols list
  So, the cols list would have [column 472, column 576, ...] where each element is a list of characters

Next, we created a sorted list which holds all the pairwise MI values of the initial columns with no bound sites (gen_mut_inf_mat)

Then, we run our algorithm:
  The while loop corresponds to each "time step"
  The first for loop takes the terms of highest MI and binds them, making a list of their indices and of their cols in the indices, cols list
    (Also removing their old individual elements)
  The next for loop removes any old MI terms from the MI list which had the newly bound columns
  Lastly, the final for loop adds the new MI terms using the bound site as an element

We continue until no MI term is >= gamma




# modules.py

This file stores some basic calculation functions that the code references.

The function shannon_entr Calculates shannon entropy of a column by adding a value to either index 0,1,2,3, or 4 of a list by virtue of the elements value and then add the end using the probabilities to convert to bits of entropy

The function joint_entr takes an unspecified amount of columns, loops through rows using zip, and generates a probability dictionary of all the occurrences to convert again to bits

The function mutual_inf takes two lists (or lists of lists) and calculates MI from joint and individual entropy

Note: the shannon_entr function is not actually in use because the method of setting up joint_entr allows for it to work even on a single column. More testing is required to determine which function is optimum for the case of a single column (I suspect the zip magic works well with many items but poorly in terms of time for just 1 item)


Thank you for reading!
