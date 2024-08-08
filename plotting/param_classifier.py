class PlotParameters:

    def __init__(self, params):

        self.object_init = {}
        self.method = {}
        self.classify(params)
    
    def classify(self, params):

        for key, value in params.items():
            if key in ["distribution", "file1", "fig_size", "fontsize"]:
                self.object_init[key] = value
            else:
                self.method[key] = value
            
