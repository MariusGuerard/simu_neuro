### Definition of the OR class containing all the characteristics about
### ORs of our experiments.


class OR:
    """
    OR is used when the flag OR_flag = True in set_parameter.py
    """

    ### Contains the affinity of all the ORs with all the odorants.
    ### (Note: It's the matrix of the Kn).

    affinities_dic = {
        "A" : {"1": 1, "2": 1, "3": 1, "4": 1, "5": 1, "6": 1},
        "B" : {"1": 1, "2": 1, "3": 1, "4": 1, "5": 1, "6": 1},
        "C" : {"1": 1, "2": 1, "3": 1, "4": 1, "5": 1, "6": 1},
        "D" : {"1": 1, "2": 1, "3": 1, "4": 1, "5": 1, "6": 1}
    }

    def __init__(self, NOr = 3, Nodor = 8, name_ORs = [], name_odors = []):


        # Contains the name of each OR.
        if len(name_ORs) > 0:
            self.name_ORs = name_ORs
        else:
            self.name_ORs = [chr(i) for i in range(ord('A'), NOr + ord('A'))]
        # Contain the name of each odor.
        if len(name_odors) > 0:
            self.name_odors = name_odors
        else:
            self.name_odors = [str(i) for i in range(Nodor)]

        # if (NOr != len(self.name_ORs)) or (Nodor != len(self.name_odors)):
        #     raise Exception("Number of ORs/odors does not fit the name")
        self.NOr = len(self.name_ORs)
        self.Nodor = len(self.name_odors)

        # Contains the affinity of all ORs for all odorants in a matrix with
        # NOr rows (for the ORs) and Nodor columns (for the odors).

        self.affinities = fill_matrix_from_lists_if_in_dict(self.name_ORs,
                                                            self.name_odors,
                                                            OR.affinities_dic)
        self.affinities_inverse = 1. / self.affinities


###########################
# HELPER FUNCTIONS FOR OR #
###########################

def fill_matrix_from_lists_if_in_dict(row_list, column_list, dic_of_dic):
    """
    For each element r_i of row_list, if dic_of_dic[r_i] exists, it look for
    each element r_j of column_list if dic_of_dic[r_i][c_j] exists. If it does
    we stock it in out_tmp[i][j].
    """
    Nrow = len(row_list)
    Ncolumn = len(column_list)
    out_tmp = []
    for i in range(Nrow):
        r_i = row_list[i]
        # print r_i
        try:
            dic_i = dic_of_dic[r_i]
            # If dic_i exists, let's continue, otherwise let's try dic_i+1
            # print dic_i
            r_i_list = []
            for j in range(Ncolumn):
                r_j = column_list[j]
                try:
                    r_i_list.append(dic_i[r_j])
                    # if dic_i[r_j] exists, let's continue, otherwise let's
                    # try dic_i[r_j+1]
                    # print dic_i[r_j]
                except:
                    r_i_list.append(np.inf)
            out_tmp.append(r_i_list)
        except:
            raise Exception("[defOR]: Characteristics of '" + str(r_i) +
                            "' are not defined.")
    return np.array(out_tmp)


#########
# DRAFT #
#########

        # aff_tmp = np.matrix((NOr, Nodor))
        # for i in range(NOr):
        #     ORi_name = self.name_ORs[i]
        #     print ORi_name
        #     try:
        #         ORi_dic = OR.affinities_dic[ORi_name]
        #         for j in range(Nodor):
        #             odorj_name = self.name_odors[j]
        #             print odorj_name
        #             try:
        #                 ORi_odorj_aff = ORi_dic[odorj_name]
        #                 aff_tmp[i][j] = ORi_odorj_aff
        #             except:
        #                 print "Odor not in OR"
        #                 pass
        #     except:
        #         print "This OR is not defined in the model  yet"
        #         pass
