import numpy as np
import pymbar

class fe_calc(object):
    """
    Computes the free energy for an APR transformation.
    """

    def __init__(self):

        self.temperature = 298.15
        self.k_b = 0.0019872041   # kcal/mol-K
        self.beta = 1/(self.k_b*self.temperature)  # Add auto updating for beta

        # List of DAT_restraints
        self.restraint_list = None

        # All raw restraint values. Organized as:
        # raw_values[phase][0:num_win][0:num_rest][np.array]
        self.raw_values = None  # Add a function to create this automatically

        # Goal is to populate these
        self.fe =     { 'attach':  None, 'pull':    None, 'release': None }
        self.fe_sem = { 'attach':  None, 'pull':    None, 'release': None }
        self.fe_matrix =     { 'attach':  None, 'pull':    None, 'release': None }
        self.fe_sem_matrix = { 'attach':  None, 'pull':    None, 'release': None }

    def _identify_active_rest(self, phase, change_param, restraint_list):
        """ Identify the restraints which are changing in the specified phase."""

        # This is a verbose approach, but it got messy when trying to be more pythonic.

        # The list of restraints which are changing
        active_rest = []
        # A boolean of which restriants in restraint_list are changing
        active_rest_bool = []

        # Iterate through the restraints and, for attach/release, check if force_constants are changing,
        # or, for pull, check if targets are changing
        for rest in restraint_list:
            if rest.phase[phase][change_param] is not None:
                # These is probably a better way to do this ... but
                # create an array in which the first element is repeated
                test_val = rest.phase[phase][change_param][0]
                test_len = rest.phase[phase][change_param].size
                test_arr = test_val*np.ones([test_len])
                # Then check if it is identical to the real array ...
                if np.allclose(rest.phase[phase][change_param], test_arr):
                    # If yes, then it's not changing, and thus not active
                    active_rest_bool.append(False)
                else:
                    # Add to active restraint list
                    active_rest.append(rest)
                    active_rest_bool.append(True)
            else:
                # If nothing is set for this phase, it's definitely not active
                active_rest_bool.append(False)

        return active_rest, active_rest_bool


    def _prepare_active_rest_data(self, phase, restraint_list, raw_values):
        """ 
        Identifies the restraints which are changing in the specified phase
        and reorders windows numerically, in case additional
        windows are added out of order after an initial run.
        """

        # Identify the restraint parameter which is expected to change for this phase
        if phase in ['attach', 'release']:
            change_param = 'force_constants'
        elif phase == 'pull':
            change_param = 'targets'
        else:
            raise Exception('Invalid phase specified: '+phase+'. Should be either attach, pull, or release')
    
        # First we identify the restraints which are changing in the specified phase.
        # active_rest is the list of active restraints, active_rest_bool is bool of
        # active restraints in the full restraint_list.
        active_rest, active_rest_bool = self._identify_active_rest(phase, change_param, restraint_list)

        # If we didn't find any active restraints, we're done
        if len(active_rest) == 0:
            return None

        # Create a window index list which matches numerical ascending order for the change_param
        reorder = np.argsort(active_rest[0].phase[phase][change_param])

        # If continuous_apr and attach/release, we need to know the order of pull too
        if active_rest[0].continuous_apr and phase in ['attach', 'release']:
            pull_active_rest, pull_active_rest_bool = self._identify_active_rest('pull', 'targets', restraint_list)
            if len(pull_active_rest) == 0:
                raise Exception('There should be active pull restraints, but none were found!')
            pull_reorder =  np.argsort(pull_active_rest[0].phase['pull']['targets'])
            pull_first_active_idx = next((i for i, rbool in enumerate(pull_active_rest_bool) if rbool is True), None)

        # Count number of windows
        num_win = len(reorder)

        # Count the number of active restraints
        num_rest = len(active_rest)

        # Identify the index of the first active restraint in the full
        # restraint_list by checking active_rest_bool. Probably don't need
        # to do this but it ensures we're checking the size of an active restraint
        first_active_idx = next((i for i, rbool in enumerate(active_rest_bool) if rbool is True), None)

        # Now count the size of data_points for each np.array in raw_values.
        # We only need to check size for first_active_idx, since they should all be same.
        # If continuous_apr and attach/release, we'll need to take the data
        # for the final attach/release window from the first/last pull window.
        if active_rest[0].continuous_apr and phase in ['attach', 'release']:
            data_points = np.zeros([num_win], np.int32)
            for win_idx in reorder[:-1]:
                data_points[win_idx] = raw_values[phase][win_idx][first_active_idx].size
            # If attach, then take data from first pull window
            if phase == 'attach':
                pull_idx = 0
            # Else, must be release, so take from last pull window
            else:
                pull_idx = -1
            data_points[-1] = raw_values['pull'][pull_reorder[pull_idx]][pull_first_active_idx].size
        # If not continuous_apr or pull, then it's more straightforward
        else:
            data_points = np.array([raw_values[phase][i][first_active_idx].size for i in reorder])

        # Get max_data_points ... will be needed for MBAR
        max_data_points = int(np.max(data_points))

        # Now create properly ordered force_constants and targets.
        force_constants = np.zeros([num_win, num_rest], np.float64)
        targets = np.zeros([num_win, num_rest], np.float64)
        for win_idx in reorder:
            for rest_idx in range(num_rest):
                force_constants[win_idx, rest_idx] = active_rest[rest_idx].phase[phase]['force_constants'][win_idx]
                if active_rest[rest_idx].mask3 is not None:
                    # Convert force constants in radians to degrees
                    force_constants[win_idx, rest_idx] *= ( (np.pi/180.0)**2 )
                targets[win_idx, rest_idx] = active_rest[rest_idx].phase[phase]['targets'][win_idx]

        # Also reorder the data to match the force_constants, targets.
        # Again, if attach/release, must take last window from pull
        ordered_values = [[] for i in range(num_win)]
        if active_rest[0].continuous_apr and phase in ['attach', 'release']:
            for win_idx in reorder[:-1]:
                for rest_idx,rest_bool in enumerate(active_rest_bool):
                    if rest_bool:
                        ordered_values[win_idx].append(raw_values[phase][win_idx][rest_idx])
            for rest_idx,rest_bool in enumerate(active_rest_bool):
                if rest_bool:
                    # Same logic as above ... must figure out which pull window to use as
                    # final attach/release window.
                    if phase == 'attach':
                        pull_idx = 0
                    else:
                        pull_idx = -1
                    ordered_values[-1].append(raw_values['pull'][pull_reorder[pull_idx]][rest_idx])
        else:
            for win_idx in reorder:
                for rest_idx,rest_bool in enumerate(active_rest_bool):
                    if rest_bool:
                        ordered_values[win_idx].append(raw_values[phase][win_idx][rest_idx])

        # Return all this in a list. Will pass to _run_mbar, which unpacks.
        # That seems inefficient, but we'll keep it for now.
        return [num_win, num_rest, data_points, max_data_points, active_rest, force_constants, targets, ordered_values]


    def _run_mbar(self, prepared_data):
        """
        Compute the free energy matrix for a series of windows. We'll follow the pymbar nomenclature for data structures.
        """

        # Unpack the prepared data
        num_win, num_rest, data_points, max_data_points, active_rest, force_constants, targets, ordered_values = prepared_data

        # Number of data points in each restraint value array
        N_k = np.array(data_points, np.int32)

        # Setup the reduced potential energy array. ie, the potential of each window's
        # coordinates in each window's potential function
        u_kln = np.zeros([num_win, num_win, max_data_points], np.float64)

        # Note, the organization of k = coordinate windows, l = potential windows
        # seems to be opposite of the documentation. But I got wrong numbers the other way around.
        for k in range(num_win):  # Coordinate windows
            for l in range(num_win): # Potential Windows
                for r,rest in enumerate(active_rest): # Restraints
                    # If this is a dihedral, we need to shift around restraint value
                    # on the periodic axis to make sure the lowest potential is used.
                    if rest.mask3 is not None and rest.mask4 is not None:
                        target = targets[l][r] # Taken from potential window, l
                        bool_list = ordered_values[k][r] < target - 180.0 # Coords from coord window, k
                        ordered_values[k][r][bool_list] += 360.0
                        bool_list = ordered_values[k][r] > target + 180.0
                        ordered_values[k][r][bool_list] -= 360.0

                # Compute the potential ... for each frame, sum the contributions for each restraint
                # Note, we multiply by kT
                u_kln[k,l,0:N_k[l]] = np.sum(self.beta*force_constants[l,:,None]*(ordered_values[k] - targets[l,:,None])**2, axis=0)

        # Setup mbar calc, and get free energies
        mbar = pymbar.MBAR(u_kln, N_k, verbose=False)
        Deltaf_ij, dDeltaf_ij, Theta_ij = mbar.getFreeEnergyDifferences(compute_uncertainty=True)

        # Add additional calc for correct uncertainties!!!

        # Put back into kcal/mol
        Deltaf_ij /= self.beta
        dDeltaf_ij /= self.beta

        return Deltaf_ij, dDeltaf_ij

    def compute_free_energy(self):
        """
        Do free energy calc.
        """

        for phase in ['attach', 'pull', 'release']: ### Need to add release, but check it's done correctly
            prepared_data = self._prepare_active_rest_data(phase, self.restraint_list, self.raw_values)
            if prepared_data:
                self.fe_matrix[phase], self.fe_sem_matrix[phase] = self._run_mbar(prepared_data)
                self.fe[phase] = self.fe_matrix[phase][0,-1]
                self.fe_sem[phase] = self.fe_sem_matrix[phase][0,-1]

                


        
        



#def compute(rest_list, rest_data, method='MBAR'):
#    
#    kB = 1.381e-23 * 6.022e23 / (4.184 * 1000.0)
#    print kB
#    temp = 298.15
#    beta = 1/(kB * temp) # beta
#
#    ### Attach Phase
#
#    #active_rest = [rest for rest in rest_list if rest.phase['attach']['force_constants'] is not None]
#    #active_rest_bool = [rest.phase['attach']['force_constants'] is not None for rest in rest_list]
#    active_rest = []
#    active_rest_bool = []
#    for rest in rest_list:
#        if rest.phase['attach']['force_constants'] is not None:
#            test_val = rest.phase['attach']['force_constants'][0]
#            test_len = rest.phase['attach']['force_constants'].size
#            test_arr = test_val*np.ones([test_len])
#            if np.allclose(rest.phase['attach']['force_constants'], test_arr):
#                active_rest_bool.append(False)
#            else:
#                active_rest.append(rest)
#                active_rest_bool.append(True)
#        else:
#            active_rest_bool.append(False)
#
#    if len(active_rest) > 0:
#        reorder = np.argsort(active_rest[0].phase['attach']['force_constants'])
#
#        num_win = active_rest[0].phase['attach']['force_constants'].size
#
#        num_rest = len(active_rest)
#
#        if active_rest[0].continuous_apr:
#            data_points = np.zeros([num_win], np.float64)
#            for i in reorder[:-1]:
#            ### Should make a check that the restraint we are sizing is actually active
#                data_points[i] = rest_data['attach'][i][0].size
#            data_points[-1] = rest_data['pull'][0][0].size
#        else:
#            data_points = np.array([rest_data['attach'][i][0].size for i in reorder])
#
#        max_data_points = int(np.max(data_points))
#
#        force_constants = np.zeros([num_win, num_rest], np.float64)
#        targets = np.zeros([num_win, num_rest], np.float64)
#        for win_idx in reorder:
#            for rest_idx in range(num_rest):
#                force_constants[win_idx, rest_idx] = active_rest[rest_idx].phase['attach']['force_constants'][win_idx]
#                if active_rest[rest_idx].mask3 is not None:
#                    force_constants[win_idx, rest_idx] *= ( (np.pi/180.0)**2 )
#                targets[win_idx, rest_idx] = active_rest[rest_idx].phase['attach']['targets'][win_idx]
#
#        rest_vals = [[] for i in range(num_win)]
#        if active_rest[0].continuous_apr:
#            for win_idx in reorder[:-1]:
#                for idx,rest_bool in enumerate(active_rest_bool):
#                    if rest_bool:
#                        rest_vals[win_idx].append(rest_data['attach'][win_idx][idx])
#            for idx,rest_bool in enumerate(active_rest_bool):
#                if rest_bool:
#                   rest_vals[-1].append(rest_data['pull'][0][idx])
#        else:
#            for win_idx in reorder:
#                for idx,rest_bool in enumerate(active_rest_bool):
#                    if rest_bool:
#                        rest_vals[win_idx].append(rest_data['attach'][win_idx][idx])
#
#        N_k = np.array(data_points, np.int32)
#
#        if method == 'MBAR':
#            u_kln = np.zeros([num_win, num_win, max_data_points], np.float64)
#        else:
#            u_kn = np.zeros([num_win], np.float64)
#
#        print rest_vals[4][0]
#
#        ### Note, the organization of k = coordinate windows, l = potential windows seems to be different than the documentation
#        for k in range(num_win):  # Coordinate windows
#            if method == 'MBAR':
#                for l in range(num_win): # Potential Windows
#                    for r,rest in enumerate(active_rest):
#                        if rest.mask3 is not None and rest.mask4 is not None:
#                            target = targets[l][r]
#                            bool_list = rest_vals[k][r] < target - 180.0
#                            rest_vals[k][r][bool_list] += 360.0
#                            bool_list = rest_vals[k][r] > target + 180.0
#                            rest_vals[k][r][bool_list] -= 360.0
#
#                    #u_kln[k,l,0:N_k[l]] += np.sum(beta*force_constants[k,:,None]*(rest_vals[l] - targets[k,:,None])**2, axis=0)
#                    u_kln[k,l,0:N_k[l]] += np.sum(beta*force_constants[l,:,None]*(rest_vals[k] - targets[l,:,None])**2, axis=0)
#                    #print k, l, np.mean(u_kln[k,l,0:N_k[l]])
#
#        mbar = pymbar.MBAR(u_kln, N_k, verbose=True)
#        Deltaf_ij, dDeltaf_ij, Theta_ij = mbar.getFreeEnergyDifferences(compute_uncertainty=True)
#
#        print "\n"
#        print Deltaf_ij
#        print "\n"
#
#        for k in range(num_win):
#            print force_constants[k][0]/force_constants[-1][0], Deltaf_ij[0][k]/beta, dDeltaf_ij[0][k]/beta
#
#
#    ### Pull Phase
#    active_rest = []
#    active_rest_bool = []
#    for rest in rest_list:
#        if rest.phase['pull']['targets'] is not None:
#            test_val = rest.phase['pull']['targets'][0]
#            test_len = rest.phase['pull']['targets'].size
#            test_arr = test_val*np.ones([test_len])
#            if np.allclose(rest.phase['pull']['targets'], test_arr):
#                active_rest_bool.append(False)
#            else:
#                active_rest.append(rest)
#                active_rest_bool.append(True)
#        else:
#            active_rest_bool.append(False)
#
#    if len(active_rest) > 0:
#        reorder = np.argsort(active_rest[0].phase['pull']['targets'])
#
#        num_win = active_rest[0].phase['pull']['targets'].size
#
#        num_rest = len(active_rest)
#
#        data_points = np.array([rest_data['pull'][i][0].size for i in reorder])
#
#        max_data_points = int(np.max(data_points))
#
#        force_constants = np.zeros([num_win, num_rest], np.float64)
#        targets = np.zeros([num_win, num_rest], np.float64)
#        for win_idx in reorder:
#            for rest_idx in range(num_rest):
#                force_constants[win_idx, rest_idx] = active_rest[rest_idx].phase['pull']['force_constants'][win_idx]
#                if active_rest[rest_idx].mask3 is not None:
#                    force_constants[win_idx, rest_idx] *= ( (np.pi/180.0)**2 )
#                targets[win_idx, rest_idx] = active_rest[rest_idx].phase['pull']['targets'][win_idx]
#
#        rest_vals = [[] for i in range(num_win)]
#        for win_idx in reorder:
#            for idx,rest_bool in enumerate(active_rest_bool):
#                if rest_bool:
#                    rest_vals[win_idx].append(rest_data['pull'][win_idx][idx])
#
#        N_k = np.array(data_points, np.int32)
#
#        if method == 'MBAR':
#            u_kln = np.zeros([num_win, num_win, max_data_points], np.float64)
#        else:
#            u_kn = np.zeros([num_win], np.float64)
#
#        ### Note, the organization of k = coordinate windows, l = potential windows seems to be different than the documentation
#        for k in range(num_win):  # Coordinate windows
#            if method == 'MBAR':
#                for l in range(num_win): # Potential Windows
#                    for r,rest in enumerate(active_rest):
#                        if rest.mask3 is not None and rest.mask4 is not None:
#                            target = targets[l][r]
#                            bool_list = rest_vals[k][r] < target - 180.0
#                            rest_vals[k][r][bool_list] += 360.0
#                            bool_list = rest_vals[k][r] > target + 180.0
#                            rest_vals[k][r][bool_list] -= 360.0
#
#                    #u_kln[k,l,0:N_k[l]] += np.sum(beta*force_constants[k,:,None]*(rest_vals[l] - targets[k,:,None])**2, axis=0)
#                    u_kln[k,l,0:N_k[l]] += np.sum(beta*force_constants[l,:,None]*(rest_vals[k] - targets[l,:,None])**2, axis=0)
#                    #print k, l, np.mean(u_kln[k,l,0:N_k[l]])
#
#        mbar = pymbar.MBAR(u_kln, N_k, verbose=True)
#        Deltaf_ij, dDeltaf_ij, Theta_ij = mbar.getFreeEnergyDifferences(compute_uncertainty=True)
#
#        for k in range(num_win):
#            print targets[k][0], Deltaf_ij[0][k]/beta, dDeltaf_ij[0][k]/beta
 