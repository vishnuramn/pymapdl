import weakref
from ansys.mapdl.core.mapdl import _MapdlCore


class Solution():
    """Collection of parameters and methods specific to the solution.

    Useful for checking the status of a solve after running
    ``mapdl.solve()`` and determining if it converged, the number of
    iterations to converge, etc.

    Also contains methods to simplify running a modal or static
    analysis.

    Examples
    --------
    Check if a solution has converged.

    >>> mapdl.solution.converged
    True

    Get the cumulative number of iterations.

    >>> mapdl.solution.n_cmit
    1.0
    """

    def __init__(self, mapdl):
        if not isinstance(mapdl, _MapdlCore):
            raise TypeError('Must be implemented from MAPDL class')
        self._mapdl_weakref = weakref.ref(mapdl)

    @property
    def _mapdl(self):
        """Return the weakly referenced instance of mapdl"""
        return self._mapdl_weakref()

    def _set_log_level(self, level):
        self._mapdl.set_log_level(level)

    @property
    def _log(self):
        return self._mapdl._log

    @property
    def time_step_size(self):
        """Time step size.

        Examples
        --------
        >>> mapdl.solution.time_step_size
        1.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'DTIME')

    @property
    def n_cmls(self):
        """Cumulative number of load steps.

        Examples
        --------
        >>> mapdl.solution.n_cmls
        1.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'NCMLS')

    @property
    def n_cmss(self):
        """Number of substeps. NOTE: Used only for static and transient analyses.

        Examples
        --------
        >>> mapdl.solution.n_cmss
        1.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'NCMSS')

    @property
    def n_eqit(self):
        """Number of equilibrium iterations.

        Examples
        --------
        >>> mapdl.solution.n_eqit
        1.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'EQIT')

    @property
    def n_cmit(self):
        """Cumulative number of iterations.

        Examples
        --------
        >>> mapdl.solution.n_cmit
        1.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'NCMIT')

    @property
    def converged(self):
        """Convergence indicator.  ``True`` when converged.

        Examples
        --------
        >>> mapdl.solution.converged
        True
        """
        return bool(self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'CNVG'))

    @property
    def mx_dof(self):
        """Maximum degree of freedom value.

        Examples
        --------
        >>> mapdl.solution.mx_dof
        -0.00020707416808476303
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'MXDVL')

    @property
    def res_frq(self):
        """Response frequency for 2nd order systems.

        Examples
        --------
        >>> mapdl.solution.res_frq
        0.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'RESFRQ')

    @property
    def res_eig(self):
        """Response eigenvalue for 1st order systems.

        Examples
        --------
        >>> mapdl.solution.res_eig
        0.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'RESEIG')

    @property
    def decent_parm(self):
        """Descent parameter.

        Examples
        --------
        >>> mapdl.solution.decent_parm
        0.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'DSPRM')

    @property
    def force_cnv(self):
        """Force convergence value.

        Examples
        --------
        >>> mapdl.solution.force_cnv
        0.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'FOCV')

    @property
    def moment_cnv(self):
        """Moment convergence value.

        Examples
        --------
        >>> mapdl.solution.moment_cnv
        0.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'MOCV')

    @property
    def heat_flow_cnv(self):
        """Heat flow convergence value.

        Examples
        --------
        >>> mapdl.solution.heat_flow_cnv
        0.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'HFCV')

    @property
    def magnetic_flux_cnv(self):
        """Magnetic flux convergence value.

        Examples
        --------
        >>> mapdl.solution.magnetic_flux_cnv
        0.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'MFCV')

    @property
    def current_segment_cnv(self):
        """Current segment convergence value.

        Examples
        --------
        >>> mapdl.solution.current_segment_cnv
        0.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'CSCV')

    @property
    def current_cnv(self):
        """Current convergence value.

        Examples
        --------
        >>> mapdl.solution.current_cnv
        0.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'CUCV')

    @property
    def fluid_flow_cnv(self):
        """Fluid flow convergence value.

        Examples
        --------
        >>> mapdl.solution.fluid_flow_cnv
        0.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'FFCV')

    @property
    def displacement_cnv(self):
        """Displacement convergence value.

        Examples
        --------
        >>> mapdl.solution.displacement_cnv
        0.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'DICV')

    @property
    def rotation_cnv(self):
        """Rotation convergence value.

        Examples
        --------
        >>> mapdl.solution.rotation_cnv
        0.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'ROCV')

    @property
    def temperature_cnv(self):
        """Temperature convergence value.

        Examples
        --------
        >>> mapdl.solution.temperature_cnv
        0.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'TECV')

    @property
    def vector_cnv(self):
        """Vector magnetic potential convergence value.

        Examples
        --------
        >>> mapdl.solution.vector_cnv
        0.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'VMCV')

    @property
    def smcv(self):
        """Scalar magnetic potential convergence value.

        Examples
        --------
        >>> mapdl.solution.smcv
        0.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'SMCV')

    @property
    def voltage_conv(self):
        """Voltage convergence value.

        Examples
        --------
        >>> mapdl.solution.voltage_conv
        0.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'VOCV')

    @property
    def pressure_conv(self):
        """Pressure convergence value.

        Examples
        --------
        >>> mapdl.solution.pressure_conv

        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'PRCV')

    @property
    def velocity_conv(self):
        """Velocity convergence value.

        Examples
        --------
        >>> mapdl.solution.velocity_conv

        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'VECV')

    @property
    def mx_creep_rat(self):
        """Maximum creep ratio.

        Examples
        --------
        >>> mapdl.solution.mx_creep_rat
        0.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'CRPRAT')

    @property
    def mx_plastic_inc(self):
        """Maximum plastic strain increment.

        Examples
        --------
        >>> mapdl.solution.mx_plastic_inc
        0.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'PSINC')

    @property
    def n_cg_iter(self):
        """Number of iterations in the PCG and symmetric JCG (non-complex version) solvers.

        Examples
        --------
        >>> mapdl.solution.n_cg_iter
        0.0
        """
        return self._mapdl.get_value('ACTIVE', 0, 'SOLU', 'CGITER')

    def run_modal_analysis(self, method='lanb', nmode='', freqb='',
                           freqe='', cpxmod='', nrmkey='', modtype='',
                           memory_option='', elcalc=False, **kwargs):
        """Run a modal analysis with basic settings analysis.

        Parameters
        ----------
        method : str
            Mode-extraction method to be used for the modal analysis.
            Defaults to lanb (block lanczos).  Must be one of the following:

            - LANB : Block Lanczos
            - LANPCG : PCG Lanczos
            - SNODE : Supernode modal solver
            - SUBSP : Subspace algorithm
            - UNSYM : Unsymmetric matrix
            - DAMP : Damped system
            - QRDAMP : Damped system using QR algorithm
            - VT : Variational Technology

        nmode : int, optional
            The number of modes to extract. The value can depend on
            the value supplied for Method. NMODE has no default and
            must be specified. If Method = LANB, LANPCG, or SNODE, the
            number of modes that can be extracted can equal the DOFs
            in the model after the application of all boundary
            conditions.

        freqb : float, optional
            The beginning, or lower end, of the frequency range of
            interest.

        freqe : float, optional
            The ending, or upper end, of the frequency range of
            interest (in Hz). The default for Method = SNODE is
            described below. The default for all other methods is to
            calculate all modes, regardless of their maximum
            frequency.

        cpxmod : str, optional
            Complex eigenmode key. Valid only when ``method='QRDAMP'``
            or ``method='unsym'``

            - AUTO : Determine automatically if the eigensolutions are
              real or complex and output them accordingly. This is
              the default for ``method='UNSYM'``.  Not supported for
              Method = QRDAMP.
            - ON or CPLX : Calculate and output complex eigenmode
              shapes.
            - OFF or REAL : Do not calculate complex eigenmode
              shapes. This is required if a mode-
              superposition analysis is intended after the
              modal analysis for Method = QRDAMP. This is the
              default for this method.

        nrmkey : bool, optional
            Mode shape normalization key.  When ``True`` (default),
            normalize the mode shapes to the mass matrix.  When False,
            Normalize the mode shapes to unity instead of to the mass
            matrix.  If a subsequent spectrum or mode-superposition
            analysis is planned, the mode shapes should be normalized
            to the mass matrix.

        modtype : str, optional
            Type of modes calculated by the eigensolver. Only
            applicable to the unsymmetric eigensolver.

            - Blank : Right eigenmodes. This value is the default.
            - BOTH : Right and left eigenmodes. The left eigenmodes are
              written to Jobname.LMODE.  This option must be
              activated if a mode-superposition analysis is intended.

        memory_option : str, optional
            Memory allocation option:

            * ``DEFAULT`` - Default Memory mode
                      Use the default memory allocation strategy for
                      the sparse solver. The default strategy attempts
                      to run in the INCORE memory mode. If there is
                      not enough available physical memory when the
                      solver starts to run in the ``INCORE`` memory
                      mode, the solver will then attempt to run in the
                      ``OUTOFCORE`` memory mode.

            * ``INCORE`` - In-core memory mode
                     Use a memory allocation strategy in the sparse
                     solver that will attempt to obtain enough memory
                     to run with the entire factorized matrix in
                     memory. This option uses the most amount of
                     memory and should avoid doing any I/O. By
                     avoiding I/O, this option achieves optimal solver
                     performance. However, a significant amount of
                     memory is required to run in this mode, and it is
                     only recommended on machines with a large amount
                     of memory. If the allocation for in-core memory
                     fails, the solver will automatically revert to
                     out-of-core memory mode.

            * ``OUTOFCORE`` - Out of core memory mode.
                        Use a memory allocation strategy in the sparse
                        solver that will attempt to allocate only
                        enough work space to factor each individual
                        frontal matrix in memory, but will store the
                        entire factorized matrix on disk. Typically,
                        this memory mode results in poor performance
                        due to the potential bottleneck caused by the
                        I/O to the various files written by the
                        solver.

        elcalc : bool, optional
            Calculate element results, reaction forces, energies, and
            the nodal degree of freedom solution.  Default ``False``.

        Returns
        -------
        response : str
            Output from MAPDL SOLVE command.

        Examples
        --------
        Modal analysis using default parameters for the first 6 modes

        >>> mapdl.modal_analysis(nmode=6)

        Notes
        -----
        For models that involve a non-symmetric element stiffness
        matrix, as in the case of a contact element with frictional
        contact, the QRDAMP eigensolver (MODOPT, QRDAMP) extracts
        modes in the modal subspace formed by the eigenmodes from the
        symmetrized eigenproblem. The QRDAMP eigensolver symmetrizes
        the element stiffness matrix on the first pass of the
        eigensolution, and in the second pass, eigenmodes are
        extracted in the modal subspace of the first eigensolution
        pass. For such non- symmetric eigenproblems, you should verify
        the eigenvalue and eigenmode results using the non-symmetric
        matrix eigensolver (MODOPT,UNSYM).

        The DAMP and QRDAMP options cannot be followed by a subsequent
        spectrum analysis. The UNSYM method supports spectrum analysis
        when eigensolutions are real.

        """
        if nrmkey:
            if nrmkey.upper() != 'OFF':
                nrmkey = 'ON'
        nrmkey = 'OFF'

        with self._mapdl.chain_commands:
            self._mapdl.slashsolu()
            self._mapdl.antype(2, 'new')
            self._mapdl.modopt(method, nmode, freqb, freqe, cpxmod, nrmkey, modtype)
            self._mapdl.bcsoption(memory_option)

            if elcalc:
                self._mapdl.mxpand(elcalc='YES')

        out = self._mapdl.solve(**kwargs)
        self._mapdl.finish(mute=True)
        return out

    def run_static_analysis(self, lg_deflect=False, **kwargs):
        """Run a new static analysis.

        Parameters
        ----------
        lg_deflect : bool, optional
            Include large-deflection effects.

        Returns
        --------
        str
            Response from MAPDL solver.

        Examples
        --------
        >>> mapdl.solution.run_static_analysis()

        Notes
        -----
        Equalivent to:

        >>> mapdl.run('/SOLU')
        >>> mapdl.antype(0)
        >>> mapdl.nlgeom(int(lg_deflect))
        >>> out = mapdl.solve()
        >>> mapdl.finish()

        """

        with self._mapdl.chain_commands:
            self._mapdl.slashsolu()
            self._mapdl.antype(0)
            self._mapdl.nlgeom(int(lg_deflect))

        out = self._mapdl.solve(**kwargs)
        self._mapdl.finish(mute=True)
        return out
