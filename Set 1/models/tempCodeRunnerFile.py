    def time_independent(self,method_fun,*args,**kwargs):
        """
        Executes a time_step ing method given a function
        Inputs:
            - method_func: Method to use function
            - stop: Stopping criteria
            - store_step: Every how many steps store data
        """
        self.data = [self.A]
        C = np.copy(self.A)
        n_count = 0
        for A_t in method_fun(C, object_ = self.object_,stop = 0.00001,*args,**kwargs):
            self.data.append(np.copy(A_t))