import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """        
    # print(f"one_gene = {one_gene}")
    # print(f"two_genes = {two_genes}")
    # print(f"have_trait = {have_trait}")
    print(f"people = {people}")

    # create data base for each person
    people_data = {}
    for person in people:        
        # find out how many genes
        if person in one_gene:
            num_gene = 1
        elif person in two_genes:
            num_gene = 2
        else:
            num_gene = 0
        
        # find out if the person has trait
        if person in have_trait:
            this_have_trait = True
        else: 
            this_have_trait = False

        # find out if the person has parents data
        if people[person]['mother'] == None:  # no parents
            this_has_parents = False
        else:
            this_has_parents = True

        people_data[person] = {"num_gene": num_gene, "trait": this_have_trait, 
                               "has_parents": this_has_parents}
    
    print(f"people_data = {people_data}")
    
    for person in people_data:        
        if people_data[person]["has_parents"] == False:
        # w/o parents data, , i.e., independent probability
            prob_gene = PROBS["gene"][people_data[person]["num_gene"]]
            prob_trait = PROBS["trait"][people_data[person]["num_gene"]][people_data[person]["trait"]]
            this_prob = prob_gene * prob_trait
            print(f"person = {person}, prob_gene = {prob_gene}, prob_trait = {prob_trait}")
        else: 
        # w/ parents data            
            match people_data[person]["num_gene"]:
                case 0:  # no gene from parents
                    # not from father                    
                    father = people[person]["father"]
                    father_gene = people_data[father]["num_gene"]
                    if father_gene == 0:
                        # father has no gene, as long as no mutation, the child won't have any
                        prob_from_father = 1 - PROBS["mutation"]
                    elif father_gene == 1:
                        # father has one gene. either from the gene with mutation or from the not-gene without mutation
                        prob_from_father = 0.5*(1 - PROBS["mutation"]) + 0.5*(PROBS["mutation"])
                    else:  # father_gene == 2
                        # father has two gene, only possibility is mutation
                        prob_from_father = PROBS["mutation"]
                    
                    # not from mother
                    mother = people[person]["mother"]
                    mother_gene = people_data[mother]["num_gene"]
                    if mother_gene == 0:
                        # mother has no gene, as long as no mutation, the child won't have any
                        prob_from_mother = 1 - PROBS["mutation"]
                    elif mother_gene == 1:
                        # mother has one gene. either from the gene with mutation or from the not-gene without mutation
                        prob_from_mother = 0.5*(1 - PROBS["mutation"]) + 0.5*(PROBS["mutation"])
                    else:  # mother_gene == 2
                        # mother has two gene, only possibility is mutation
                        prob_from_mother = PROBS["mutation"]

                    prob_gene = prob_from_father * prob_from_mother
                    prob_trait = PROBS["trait"][people_data[person]["num_gene"]][people_data[person]["trait"]]
                    this_prob = prob_gene * prob_trait
                    print(f"person = {person}, prob_gene = {prob_gene}, prob_trait = {prob_trait}")
                # TO-DO
                # case 1:
                # case 2:

    


    raise NotImplementedError


def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    raise NotImplementedError


def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    raise NotImplementedError


if __name__ == "__main__":
    main()
