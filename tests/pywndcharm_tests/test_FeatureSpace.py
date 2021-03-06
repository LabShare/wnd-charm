"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                               
 Copyright (C) 2015 National Institutes of Health 

    This library is free software; you can redistribute it and/or              
    modify it under the terms of the GNU Lesser General Public                 
    License as published by the Free Software Foundation; either               
    version 2.1 of the License, or (at your option) any later version.         
                                                                               
    This library is distributed in the hope that it will be useful,            
    but WITHOUT ANY WARRANTY; without even the implied warranty of             
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          
    Lesser General Public License for more details.                            
                                                                               
    You should have received a copy of the GNU Lesser General Public           
    License along with this library; if not, write to the Free Software        
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA  
                                                                               
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                               
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Written by:  Christopher Coletta (github.com/colettace)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""

import sys
if sys.version_info < (2, 7):
    import unittest2 as unittest
else:
    import unittest


#for i, _ in enumerate(sys.path):
#    if '/usr/local/lib/python2.7/site-packages/wndcharm-0.9.13' in _:
#        del sys.path[i]
#        print "deleting", _
#    else:
#        print i, _
#from wndcharm import diagnostics
#print diagnostics

import re
import numpy as np

from os.path import sep, dirname, realpath, join
from tempfile import mkdtemp
from shutil import rmtree

pychrm_test_dir = dirname( realpath( __file__ ) ) #WNDCHARM_HOME/tests/pychrm_tests
wndchrm_test_dir = dirname( pychrm_test_dir ) + sep + 'wndchrm_tests'

from wndcharm.FeatureSpace import FeatureSpace
from wndcharm.FeatureWeights import FisherFeatureWeights, PearsonFeatureWeights
from wndcharm.FeatureSpacePrediction import FeatureSpaceClassification, \
       FeatureSpaceRegression
from wndcharm.utils import compare
from wndcharm.FeatureVector import FeatureVector

class TestFeatureSpace( unittest.TestCase ):
    """
    The FeatureSet is the workhorse object in WND-CHARM.
    """
    maxDiff = None
    test_fit_path = join( wndchrm_test_dir,'test-l.fit' )
    test_normalized_fit_path = join( wndchrm_test_dir, 'test_fit-l-normalized.fit' )

    # --------------------------------------------------------------------------
    #@unittest.skip('')
    def test_ContinuousFitOnFit( self ):
        from wndcharm.ArtificialFeatureSpace import CreateArtificialFeatureSpace_Discrete

        fs_discrete = CreateArtificialFeatureSpace_Discrete( n_samples=1000, n_classes=10,
                num_features_per_signal_type=30, noise_gradient=5, initial_noise_sigma=10,
                n_samples_per_group=1, interpolatable=True)

        tempdir = mkdtemp()
        path_to_fit = tempdir + sep + 'Artificial.fit'
 
        try:
          fs_discrete.ToFitFile( path_to_fit )
          fs_continuous = FeatureSpace.NewFromFitFile(
                  path_to_fit, discrete=False )

          fs_continuous.Normalize( quiet=True )
          fw_reduced = PearsonFeatureWeights.NewFromFeatureSpace( fs_continuous ).Threshold()
          fs_reduced = fs_continuous.FeatureReduce( fw_reduced )
          batch_result = FeatureSpaceRegression.NewMultivariateLinear( 
                  fs_reduced, fw_reduced, quiet=True )

        finally:
          rmtree( tempdir )

    # --------------------------------------------------------------------------
    #@unittest.skip('')
    def test_DiscreteTrainTestSplitNoTiling( self ):
        """Uses binucleate test set"""

        fitfilepath = wndchrm_test_dir + sep + 'test-l.fit'
        fs = FeatureSpace.NewFromFitFile( fitfilepath )

        from numpy.random import RandomState
        prng = RandomState(42)
        full_train, full_test = fs.Split( random_state=prng, quiet=True )
        full_train.Normalize( quiet=True )
        reduced_fw = FisherFeatureWeights.NewFromFeatureSpace( full_train ).Threshold()
        reduced_train = full_train.FeatureReduce( reduced_fw )

        reduced_test = full_test.FeatureReduce( reduced_fw )
        reduced_test.Normalize( reduced_train, quiet=True )

        batch_result = FeatureSpaceClassification.NewWND5( reduced_train,
            reduced_test, reduced_fw, quiet=True )


    # --------------------------------------------------------------------------
    #@unittest.skip('')
    def test_DiscreteTrainTestSplitWithTiling( self ):
        """Uses a curated subset of the IICBU 2008 Lymphoma dataset, preprocessed as follows:
        auto-deconvolved, eosin channel only, tiled 5x6, 3 classes, 10 imgs per class,
        300 samples per class.
        """

        # Inflate the zipped test fit into a temp file
        import zipfile
        zipped_file_path = pychrm_test_dir + sep + 'lymphoma_iicbu2008_subset_EOSIN_ONLY_t5x6_v3.2features.fit.zip'
        zf = zipfile.ZipFile( zipped_file_path, mode='r' )
        tempdir = mkdtemp()
        zf.extractall( tempdir )

        try:
            fitfilepath = tempdir + sep + zf.namelist()[0]
            #fs = FeatureSet.NewFromFitFile( fitfilepath  )
            fs = FeatureSpace.NewFromFitFile( fitfilepath, tile_num_rows=5, tile_num_cols=6 )
            from numpy.random import RandomState
            prng = RandomState(42)
            train, test = fs.Split( random_state=prng, quiet=True )
            train.Normalize( inplace=True, quiet=True )
            fw = FisherFeatureWeights.NewFromFeatureSpace( train ).Threshold()
            train.FeatureReduce( fw, inplace=True )
            test.FeatureReduce( fw, inplace=True ).Normalize( train, inplace=True, quiet=True )

        finally:
            rmtree( tempdir )

    # --------------------------------------------------------------------------
    #@unittest.skip('')
    def test_ContinuousTrainTestSplitWithTiling( self ):
        """Uses a synthetic preprocessed as follows: 500 total samples, 25 tiles per group
        240 total features"""

        from wndcharm.ArtificialFeatureSpace import CreateArtificialFeatureSpace_Continuous

        fs = CreateArtificialFeatureSpace_Continuous( n_samples=500,
                num_features_per_signal_type=20, n_samples_per_group=25 )

        from numpy.random import RandomState
        prng = RandomState(42)
        #fs.Print( verbose=True )
        #print "\n\n\n********************\n\n\n"
        train, test = fs.Split( random_state=prng, quiet=True )
        #full_train.Print( verbose=True )
        #full_test.Print( verbose=True )

        train.Normalize( inplace=True, quiet=True )
        fw = PearsonFeatureWeights.NewFromFeatureSpace( train ).Threshold()
        train.FeatureReduce( fw, inplace=True )
        test.FeatureReduce( fw, inplace=True ).Normalize( train, inplace=True, quiet=True )

    # --------------------------------------------------------------------------
    #@unittest.skip('')
    def test_SplitOptions( self ):
        from wndcharm.ArtificialFeatureSpace import CreateArtificialFeatureSpace_Discrete

        fs_discrete = CreateArtificialFeatureSpace_Discrete( n_samples=1000, n_classes=10,
                num_features_per_signal_type=30, noise_gradient=5, initial_noise_sigma=10,
                n_samples_per_group=1, interpolatable=True, random_state=42)

        # default
        train_set, test_set = fs_discrete.Split( random_state=42, quiet=True )
        self.assertEqual( train_set.shape, (750, 600) )
        self.assertEqual( test_set.shape, (250, 600) )

        # Supposed to only return single FeatureSpace instead of 2-tuple of FeatureSpace
        # when setting test_size = 0
        i = 50
        retval = fs_discrete.Split( train_size=i, test_size=0, random_state=42, quiet=True )
        self.assertEqual( type(retval), FeatureSpace )
        self.assertEqual( retval.num_samples, i * fs_discrete.num_classes )

        # dummyproofing

        self.assertRaises( ValueError, fs_discrete.Split, train_size='trash' )
        self.assertRaises( ValueError, fs_discrete.Split, train_size=1.1 )
        self.assertRaises( ValueError, fs_discrete.Split, test_size='trash' )
        self.assertRaises( ValueError, fs_discrete.Split, test_size=1.1 )

        # What if the feature set number of groups within a class are less than called for
        # when specifying by integer?
        self.assertRaises( ValueError, test_set.Split, test_size=25 )

        # What happens when input fs has unbalanced classes, some of which have enough
        # to satisfy train_size/test_size params, and some don't
        remove_these = range(250,300) + range(700,750)
        fs_class_2_and_7_smaller = \
              fs_discrete.SampleReduce( leave_out_sample_group_ids=remove_these )

        self.assertRaises( ValueError, fs_class_2_and_7_smaller.Split, train_size=80,
                           test_size=20 )

        # Test balanced_classes:
        train_fs, test_fs = fs_class_2_and_7_smaller.Split()
        # Training set number rounds down (apparently).
        from math import floor
        expected_num_samps_per_train_class = int( floor(50*0.75) )
        expected_num_samps_per_test_class = 50 - expected_num_samps_per_train_class

        err_msg = "Balanced classes {} set split error, class {}, expected {}, got {}"
        for i, (n_train, n_test) in enumerate( zip( train_fs.class_sizes, test_fs.class_sizes )):
            self.assertEqual( n_train, expected_num_samps_per_train_class, msg=\
                    err_msg.format( "TRAIN", i, expected_num_samps_per_train_class, n_train  ) )
            self.assertEqual( n_test, expected_num_samps_per_test_class, msg=\
                    err_msg.format( "TEST", i, expected_num_samps_per_test_class, n_test ) )

    # --------------------------------------------------------------------------
    #@unittest.skip('')
    def test_SampleReduce( self ):
        from wndcharm.ArtificialFeatureSpace import CreateArtificialFeatureSpace_Discrete

        n_classes = 10
        #========================================================
        # Section 1: LEAVE IN, Untiled Discrete (w/ classes) FeatureSpace instances
        fs_discrete = CreateArtificialFeatureSpace_Discrete( n_samples=1000, n_classes=n_classes,
                num_features_per_signal_type=30, noise_gradient=5, initial_noise_sigma=10,
                n_samples_per_group=1, interpolatable=True)

        # Reduce to 9 classes from 10, one sample per class
        # Drop the last class:
        desired = range(50, 950, 100)

        A = fs_discrete.SampleReduce( desired )
        # Further reduce to 8 classes
        A.RemoveClass( "FakeClass-055.6", inplace=True )

        correct_samplenames = ['FakeClass-100.0_050', 'FakeClass-077.8_050', 'FakeClass-033.3_050', 'FakeClass-011.1_050', 'FakeClass+011.1_050', 'FakeClass+033.3_050', 'FakeClass+055.6_050', 'FakeClass+077.8_050']
        #The actual alphanumeric sort order is different from the value sort order
        #correct_samplenames = ['FakeClass+011.1_050', 'FakeClass+033.3_050', 'FakeClass+055.6_050', 'FakeClass+077.8_050', 'FakeClass-011.1_050', 'FakeClass-033.3_050', 'FakeClass-077.8_050', 'FakeClass-100.0_050']
        self.assertEqual( correct_samplenames, A._contiguous_sample_names )

        correct_classnames = ['FakeClass-100.0', 'FakeClass-077.8', 'FakeClass-033.3', 'FakeClass-011.1',  'FakeClass+011.1', 'FakeClass+033.3', 'FakeClass+055.6', 'FakeClass+077.8'  ]
        #correct_classnames = ['FakeClass+011.1', 'FakeClass+033.3', 'FakeClass+055.6', 'FakeClass+077.8', 'FakeClass-011.1', 'FakeClass-033.3', 'FakeClass-077.8', 'FakeClass-100.0']
        self.assertEqual( correct_classnames, A.class_names )
        del A

        #========================================================
        # Section 2: LEAVE OUT, UNTiled Feature sets, Discrete FeatureSpace instances

        UNdesired = range(50, 950, 100)
        C = fs_discrete.SampleReduce( leave_out_sample_group_ids=UNdesired )
        self.assertEqual( C.num_samples, fs_discrete.num_samples - len( UNdesired ) )

        # Single integers for leave_out_list is ok
        UNdesired = 50
        C = fs_discrete.SampleReduce( leave_out_sample_group_ids=UNdesired )
        self.assertEqual( C.num_samples, fs_discrete.num_samples - 1 )
        del C

        #========================================================
        # Section 3: LEAVE IN, Tiled Feature sets, Discrete FeatureSpace instances
        num_tiles = 4
        fs_discrete = CreateArtificialFeatureSpace_Discrete( n_samples=1000, n_classes=n_classes,
                num_features_per_signal_type=30, noise_gradient=5, initial_noise_sigma=10,
                n_samples_per_group=num_tiles, interpolatable=True)

        desired = range(5, 95, 10) # Rearrange into 9 classes
        D = fs_discrete.SampleReduce( desired )
        # Total num samples should be 9 classes, 1 sample group per class, 4 tiles per SG = 36
        self.assertEqual( num_tiles * len( desired ), D.num_samples )
        del D

        #========================================================
        # Section 4: LEAVE OUT, WITH Tiled Feature sets, Discrete FeatureSpace instances

        # You can't leave out a sample group that doesn't exist
        UNdesired = range(50000, 50010)
        self.assertRaises( ValueError, fs_discrete.SampleReduce,
                leave_out_sample_group_ids=UNdesired )

        # Can't leave out trash
        UNdesired = ['foo', 'bar']
        self.assertRaises( TypeError, fs_discrete.SampleReduce,
                leave_out_sample_group_ids=UNdesired )

        # This input is ok:
        UNdesired = range(5, 95, 10)
        E = fs_discrete.SampleReduce( leave_out_sample_group_ids=UNdesired )
        self.assertEqual( E.num_samples, fs_discrete.num_samples - len( UNdesired ) * num_tiles )
        del E

        #========================================================
        # Section 5: LEAVE IN, Untiled Continuous FeatureSpace instances
        from wndcharm.ArtificialFeatureSpace import CreateArtificialFeatureSpace_Continuous

        fs_cont = CreateArtificialFeatureSpace_Continuous( n_samples=1000,
                num_features_per_signal_type=30, noise_gradient=5, initial_noise_sigma=10,
                n_samples_per_group=1)

        # dummyproof
        desired = ['foo', 'bar']
        self.assertRaises( TypeError, fs_cont.SampleReduce, desired )

        desired = range(50, 950)
        F = fs_cont.SampleReduce( desired )
        self.assertEqual( F.num_samples, len(desired) )
        del F

        #========================================================
        # Section 6: LEAVE OUT, Untiled Continuous FeatureSpace instances

        UNdesired = range(50, 950)
        G = fs_cont.SampleReduce( leave_out_sample_group_ids=UNdesired )
        self.assertEqual( G.num_samples, fs_cont.num_samples - len(UNdesired) )
        del G

        # single int is ok
        H = fs_cont.SampleReduce( leave_out_sample_group_ids=998 )
        self.assertEqual( H.num_samples, fs_cont.num_samples - 1 )
        del H

        #========================================================
        # Section 7: LEAVE IN, TILED Continuous FeatureSpace instances

        fs_cont = CreateArtificialFeatureSpace_Continuous( n_samples=1000,
                num_features_per_signal_type=30, noise_gradient=5, initial_noise_sigma=10,
                n_samples_per_group=num_tiles)

        desired = range(50, 95)
        I = fs_cont.SampleReduce( desired )
        self.assertEqual( I.num_samples, len(desired) * num_tiles )
        del I

        # single int is ok, ALTHOUGH NOT SURE WHY YOU'D EVER WANT A FS WITH A SINGLE SAMPLE
        J = fs_cont.SampleReduce( 98 )
        self.assertEqual( J.num_samples, num_tiles )
        del J

        #========================================================
        # Section 8: LEAVE OUT, TILED Continuous FeatureSpace instances

        UNdesired = range(50, 95)
        K = fs_cont.SampleReduce( leave_out_sample_group_ids=UNdesired )
        self.assertEqual( K.num_samples, fs_cont.num_samples - len(UNdesired) * num_tiles )
        del K

        # single int is ok
        L = fs_cont.SampleReduce( leave_out_sample_group_ids=98 )
        self.assertEqual( L.num_samples, fs_cont.num_samples - num_tiles  )
        del L

    # --------------------------------------------------------------------------
    #@unittest.skip('')
    def test_RemoveClass( self ):
        from wndcharm.ArtificialFeatureSpace import CreateArtificialFeatureSpace_Discrete

        n_classes = 10
        fs = CreateArtificialFeatureSpace_Discrete( n_samples=1000, n_classes=n_classes,
                num_features_per_signal_type=30, noise_gradient=5, initial_noise_sigma=10,
                n_samples_per_group=1, interpolatable=True )
        orig_num_features = fs.num_features

        #print "BEFORE:", str( fs )
        fs.RemoveClass( 5, inplace=True )
        #print "AFTER:", str( fs )

        self.assertEqual( fs.num_classes, n_classes - 1 )

        # Now, the list of classes are:
        # ['FakeClass-100', 'FakeClass-77.78', 'FakeClass-55.56', 'FakeClass-33.33', 'FakeClass-11.11', 'FakeClass33.33', 'FakeClass55.56', 'FakeClass77.78', 'FakeClass100']
        # "FakeClass33.33" is the 5th class, so essentially we're trying to remove the
        # 5th class twice:

        fs.RemoveClass( "FakeClass+033.3" , inplace=True )
        #print "AFTER AFTER:", str( fs )

        self.assertEqual( fs.num_classes, n_classes - 2 )
        # Each class has 100 samples, and we've removed 2 classes:
        self.assertEqual( fs.num_samples, 800 )
        self.assertEqual( fs.shape[0], 800 )
        self.assertEqual( fs.num_features, orig_num_features )

        self.assertRaises( ValueError, fs.RemoveClass, class_token='trash' )
        self.assertRaises( IndexError, fs.RemoveClass, class_token=10 )

    # --------------------------------------------------------------------------
    #@unittest.skip('')
    def test_SamplesUnion( self ):
        from wndcharm.ArtificialFeatureSpace import CreateArtificialFeatureSpace_Discrete

        n_classes = 2
        fs1 = CreateArtificialFeatureSpace_Discrete( n_samples=20, n_classes=n_classes,
                num_features_per_signal_type=30, noise_gradient=5, initial_noise_sigma=10,
                n_samples_per_group=1, interpolatable=True )

        fitfile_path = wndchrm_test_dir + sep + 'test-l.fit'
        fs2 = FeatureSpace.NewFromFitFile( fitfile_path )

        self.assertRaises( ValueError, fs1.SamplesUnion, other_fs=fs2 )

        fs3 = CreateArtificialFeatureSpace_Discrete( n_samples=20, n_classes=n_classes,
                num_features_per_signal_type=30, noise_gradient=5, initial_noise_sigma=10,
                n_samples_per_group=1, interpolatable=True )

        joined_fs = fs1 + fs3

        self.assertEqual( n_classes, joined_fs.num_classes )

    # --------------------------------------------------------------------------
    #@unittest.skip('')
    def test_NewFromFileOfFiles( self ):
        """Pulls in the lymphoma eosin histology 5x6 tiled featureset via .sig files."""

        # Types of files containing features:
        # FIT: contains an entire FeatureSpace definition including features.
        # FOF: "File Of Files" containing a FeatureSpace structure definition only,
        #      listing paths to files of pre-calculated features (.sig) or the
        #      tiff images themselves so features can be calculated
        # SIG: A text file containing pre-calculated features for a single sample.

        # Test dataset: subset of the IICBU2008 lymphoma dataset. 2 channels (H+E),
        #    3 classes ('CLL', 'FL', 'MCL'), 10 images per class per channel,
        #    5x6 tiling grid = 30 samples per image resulting in 
        #    2 x 3 x 10 X 30 = 1800 total samples available

        # Files containing features included in this test suite:
        # 1. lymphoma_iicbu2008_subset_EOSIN_ONLY_t5x6_v3.2features.fit.zip:
        #    A zip archive containing a single FIT file with features pre-calculated.
        # 2. lymphoma_iicbu2008_subset_HE_t5x6_v3.2features_SIGFILES.zip:
        #    Contains 1800 SIG files, plus 4 FOF files (items 2-5 below):
        #       "lymphoma_iicbu2008_subset_EOSIN_ONLY_images.fof.tsv"
        #       "lymphoma_iicbu2008_subset_EOSIN_ONLY_sigfiles_t5x6-l.fof.tsv"
        #       "lymphoma_iicbu2008_subset_2CHAN_HE_images.fof.tsv"
        #       "lymphoma_iicbu2008_subset_2CHAN_HE_sigfiles_t5x6-l.fof.tsv"

        # List of possible feature sources:
        #    1. Single channel FIT (Eosin only)
        #    2. Single channel FOF (Eosin only) referencing to 30 tiffs (requires global sampling options -t5x6 -l to grab sigs)
        #    3. Single channel FOF (Eosin only) referencing 900 sig files
        #    4. Double channel FOF (Eosin+Haemotoxylin) referencing 60 tiffs (requires global sampling options -t5x6 -l to grab sigs)
        #    5. Double channel FOF (Eosin+Haemotoxylin) referencing 1800 sig files.

        #=============================================
        # BEGIN CODE TO CREATE TESTDATA ZIP PACKAGE

        #import zipfile
        #import zlib
        #path = '/Users/chris/src/wnd-charm/tests/pywndcharm_tests/TESTDATA_lymphoma_iicbu2008_subset_HE_t5x6_v3.2features_SIGFILES.zip'
        #zf = zipfile.ZipFile( path, mode='w' )
        #import os
        #classes = 'CLL', 'FL', 'MCL',
        #channels = 'haemotoxylin', 'eosin'
        #from collections import defaultdict
        #sig_tracker = defaultdict(int)
        #samplegroupid_tracker = {}
        #samplegroup_counter = 0
        #
        #eosin_tif_fof = [] # 30 lines
        #eosin_sig_fof = [] # 900 lines
        #double_tif_fof = [] # 30 lines, 2 feature set columns
        #double_sig_fof = [] # 900 lines, 2 feature set columns
        #
        #for _channel in channels:
        #    zf.write( './' + _channel, compress_type=zipfile.ZIP_DEFLATED )
        #    for _class in classes:
        #        zf.write( './' + _channel + '/' + _class, compress_type=zipfile.ZIP_DEFLATED )
        #        for root, dirs, files in os.walk( _channel + '/' + _class ):
        #            for _file in files:
        #                if _file.endswith( '.tif' ):
        #                    # Strip off the _H.tif or _E.tif
        #                    samplename = _file[:-6]
        #                    eosinpath = './eosin/' + _class + '/' + samplename + '_E.tif'
        #                    haemopath = './haemotoxylin/' + _class + '/' + samplename + '_H.tif'
        #                    if _channel == 'eosin':
        #                        eosin_tif_fof.append( eosinpath + '\t' + _class )
        #                        double_tif_fof.append( samplename + '\t' + _class + '\t' + eosinpath + '\t{\tchannel\t=\teosin\t}\t' + haemopath + '\t{\tchannel\t=\thaemotoxylin\t}')
        #                elif _file.endswith( '.sig' ):
        #                    zf.write( './' + _channel + '/' + _class + '/' + _file, compress_type=zipfile.ZIP_DEFLATED )
        #                    if _channel == 'eosin':
        #                        # Strip off the _H-t5x6_0_0-l.sig
        #                        samplename = _file[:-17] + '.tif'
        #                        eosinpath = './eosin/' + _class + '/' + _file
        #                        haemopath = './haemotoxylin/' + _class + '/' + _file.replace( '_E-t5x6_', '_H-t5x6_' )
        #                        # count samples from 0:
        #                        samplesequenceid = str( sig_tracker[ samplename ] )
        #                        sig_tracker[ samplename ] += 1
        #                        if samplename not in samplegroupid_tracker:
        #                            samplegroupid_tracker[ samplename ] = samplegroup_counter
        #                            samplegroup_counter += 1
        #                        samplegroupid = str( samplegroupid_tracker[ samplename ] )
        #                        eosin_sig_fof.append( eosinpath + '\t' + _class )
        #                        double_sig_fof.append( samplename + '\t' + _class + '\t' + eosinpath + '\t{\tchannel\t=\teosin\t;\tsamplegroupid\t=\t' + samplegroupid + '\t;\tsamplesequenceid\t=\t' + samplesequenceid + '\t}\t' + haemopath + '\t{\tchannel\t=\thaemotoxylin\t;\tsamplegroupid\t=\t' + samplegroupid + '\t;\tsamplesequenceid\t=\t' + samplesequenceid + '\t}\t')
        #
        #fof_dir = '/Users/chris/src/wnd-charm/tests/pywndcharm_tests/'
        #with open( 'lymphoma_iicbu2008_subset_EOSIN_ONLY_images.fof.tsv', 'w') as out:
        #    for _ in eosin_tif_fof:
        #        out.write( _ + '\n')
        #with open( 'lymphoma_iicbu2008_subset_EOSIN_ONLY_sigfiles_t5x6-l.fof.tsv', 'w') as out:
        #    for _ in eosin_sig_fof:
        #        out.write( _ + '\n')
        #with open( 'lymphoma_iicbu2008_subset_2CHAN_HE_images.fof.tsv', 'w') as out:
        #    for _ in double_tif_fof:
        #        out.write( _ + '\n')
        #with open( 'lymphoma_iicbu2008_subset_2CHAN_HE_sigfiles_t5x6-l.fof.tsv', 'w') as out:
        #    for _ in double_sig_fof:
        #        out.write( _ + '\n')
        #zf.write( './' + 'lymphoma_iicbu2008_subset_EOSIN_ONLY_images.fof.tsv', compress_type=zipfile.ZIP_DEFLATED )
        #zf.write( './' + 'lymphoma_iicbu2008_subset_EOSIN_ONLY_sigfiles_t5x6-l.fof.tsv', compress_type=zipfile.ZIP_DEFLATED )
        #zf.write( './' + 'lymphoma_iicbu2008_subset_2CHAN_HE_images.fof.tsv', compress_type=zipfile.ZIP_DEFLATED )
        #zf.write( './' + 'lymphoma_iicbu2008_subset_2CHAN_HE_sigfiles_t5x6-l.fof.tsv', compress_type=zipfile.ZIP_DEFLATED )
        #zf.printdir()
        #zf.close()

        # END CODE TO CREATE TESTDATA ZIP PACKAGE
        #=============================================

        # Inflate the zipped test fit into a temp file
        import zipfile
        
        zipped_file_path = pychrm_test_dir + sep + 'lymphoma_iicbu2008_subset_HE_t5x6_v3.2features_SIGFILES.zip'
        zf1 = zipfile.ZipFile( zipped_file_path, mode='r' )
        tempdir = mkdtemp()
        zf1.extractall( tempdir )

        # for comparison:
        zf2 = zipfile.ZipFile( pychrm_test_dir + sep + 'lymphoma_iicbu2008_subset_EOSIN_ONLY_t5x6_v3.2features.fit.zip', mode='r')
        zf2.extractall( tempdir )

        try:
            kwargs = {}
            kwargs['pathname'] = tempdir + sep + 'lymphoma_iicbu2008_subset_EOSIN_ONLY_sigfiles_t5x6-l.fof.tsv'
            kwargs['quiet'] = True
            # sampling opts: -l -t5x6 implies 5 columns and 6 rows ... I know it's weird.
            kwargs['long'] = True
            kwargs['tile_num_rows'] = 6
            kwargs['tile_num_cols'] = 5
            fs_fof = FeatureSpace.NewFromFileOfFiles( **kwargs )

            kwargs['pathname'] = tempdir + sep + 'lymphoma_iicbu2008_subset_eosin_t5x6_v3.2features.fit'
            fs_fit = FeatureSpace.NewFromFitFile( **kwargs )

            # Fit file has less significant figures than Signature files, and it's not
            # consistent how many there are. Seems like fit file just lops off numbers
            # at the end. Example: (signatures on top, fit on bottom)
            #
            # Example:
            # -  17.232246,  # sig
            # ?         --
            #
            # +  17.2322,    # fit
            # -  -63.549056, # sig
            # ?         ^^^
            #
            # +  -63.5491,   # fit
            # ?         ^
            #
            # -  223.786977, # sig
            # ?        ---
            #
            # +  223.787,    # fit

            # More of the same:
            #(Pdb) fs_fof.data_matrix[0,-5:]
            #array([   0.935442,   14.005003,  -43.562076,  127.394914,    0.628772])
            #(Pdb) fs_fit.data_matrix[0,-5:]
            #array([   0.935442,   14.005   ,  -43.5621  ,  127.395   ,    0.628772])

            # default is rtol=1e-07, atol=0
            #np.testing.assert_allclose( actual=fs_fit.data_matrix, desired=fs_fof.data_matrix,
            #        rtol=1e-03, atol=0 )
            #np.testing.assert_array_almost_equal_nulp( fs_fit.data_matrix, fs_fof.data_matrix )
            for row_num, (fit_row, fof_row) in enumerate( zip( fs_fit.data_matrix, fs_fof.data_matrix )):
                retval = compare( fit_row, fof_row )
                if retval == False:
                    print "error in sample row", row_num
                    print "FIT: ", fs_fit._contiguous_sample_names[row_num], "FOF", fs_fof._contiguous_sample_names[row_num]
                self.assertTrue( retval )


            # Test sorting; scramble the FOF then load and check:

            sorted_fof = tempdir + sep + \
                    'lymphoma_iicbu2008_subset_EOSIN_ONLY_sigfiles_t5x6-l.fof.tsv'

            with open( sorted_fof) as fof:
                lines = fof.readlines()

            from random import shuffle
            shuffle(lines)

            unsorted_fof = tempdir + sep + \
                    'lymphoma_iicbu2008_subset_EOSIN_ONLY_sigfiles_t5x6-l_UNSORTED.fof.tsv'

            with open( unsorted_fof, 'w' ) as fof:
                for line in lines:
                    fof.write( line )

            kwargs = {}
            kwargs['pathname'] = unsorted_fof
            kwargs['quiet'] = True
            # sampling opts: -l -t5x6 implies 5 columns and 6 rows ... I know it's weird.
            kwargs['long'] = True
            kwargs['tile_num_rows'] = 6
            kwargs['tile_num_cols'] = 5
            fs_fof = FeatureSpace.NewFromFileOfFiles( **kwargs )
            # Check again
            for row_num, (fit_row, fof_row) in enumerate( zip( fs_fit.data_matrix, fs_fof.data_matrix )):
                retval = compare( fit_row, fof_row )
                if retval == False:
                    print "error in sample row", row_num
                    print "FIT: ", fs_fit._contiguous_sample_names[row_num], "FOF", fs_fof._contiguous_sample_names[row_num]
                self.assertTrue( retval )

            # TESTING TAKE TILES:
            self.assertRaises( ValueError, fs_fof.TakeTiles, tuple() )
            self.assertRaises( ValueError, fs_fof.TakeTiles, (45, 46, 47,) )
            self.assertRaises( TypeError, fs_fof.TakeTiles, 'crap' )

            # take middle 4
            wanted_tiles = ( 14, 15, 20, 21 )

            took = fs_fof.TakeTiles( wanted_tiles, inplace=False )
            num_sample_groups = len( set( fs_fof._contiguous_sample_group_ids ) )
            self.assertEqual( took.num_samples_per_group, len( wanted_tiles ) )
            self.assertEqual( took.num_samples, len( wanted_tiles ) * num_sample_groups )

#            mid4 = 'lymphoma_iicbu2008_subset_EOSIN_ONLY_sigfiles_MIDDLE_4_TILES_t5x6-l.fof.tsv'
#            # fake out wndcharm by putting empty tiffs in the temp dir
#            # we don't need them, the sigs are in there already.
#            with open( mid4) as fof:
#                lines = fof.readlines()
#                names, classes, paths, opts = zip( *[ _.split('\t') for _ in lines ] )
#                for _path in paths:
#                    with open( tempdir + sep + _path, 'w' ):
#                        pass
#            took_via_fof = FeatureSpace.NewFromFileOfFiles( mid4, num_samples_per_group=4 )
#
#            for row_num, (fit_row, fof_row) in enumerate( zip( took.data_matrix, took_via_fof.data_matrix )):
#                retval = compare( fit_row, fof_row )
#                if retval == False:
#                    print "error in sample row", row_num
#                    print "FIT: ", took._contiguous_sample_names[row_num], "FOF", took_via_fof._contiguous_sample_names[row_num]
#                self.assertTrue( retval )


        finally:
            rmtree( tempdir )

    # --------------------------------------------------------------------------
    @unittest.skip('')
    def test_NewFromFileOfFiles_PROFILE( self ):
        """Profile the loading of a TSV"""
        import cProfile as pr
        #import profile as pr
        import tempfile
        import pstats
        prof = tempfile.NamedTemporaryFile()
        cmd = 'fs = FeatureSpace.NewFromFileOfFiles( "/home/maltsevav/Documents/Notebooks/10-scale-12-tiles-grayscale-50.tsv", long=True, num_samples_per_group=12 )' 
        try:
            pr.runctx( cmd, globals(), locals(), prof.name)
        except StandardError as e:
            print e
        p = pstats.Stats(prof.name)
        #fs = FeatureSpace.NewFromFileOfFiles( "/home/maltsevav/Documents/Notebooks/10-scale-12-tiles-grayscale-50.tsv", long=True, num_samples_per_group=12 )
        p.sort_stats('time').print_stats(50)
        print '\n\n\n**********************************\n\n\n'
        p.sort_stats('cumtime').print_stats(50)
        prof.close()
    # --------------------------------------------------------------------------
    #@unittest.skip('')
    def test_NewFromDirectory( self ):
        """"""

        ref_sig_path = 'lymphoma_eosin_channel_MCL_test_img_sj-05-3362-R2_001_E-t6x5_5_4-l.sig'
        ref_fv = FeatureVector.NewFromSigFile( pychrm_test_dir + sep + ref_sig_path )
        from shutil import copy
        tempdir = mkdtemp()
        img_filename = "lymphoma_eosin_channel_MCL_test_img_sj-05-3362-R2_001_E-t6x5_5_4-l.tiff"
        orig_img_filepath = pychrm_test_dir + sep + img_filename
        copy( orig_img_filepath, tempdir )
        try:
            fs = FeatureSpace.NewFromDirectory( tempdir, quiet=False )
            self.assertTrue( compare( fs.data_matrix[0], ref_fv.values ) )
            #from numpy.testing import assert_allclose
            #assert_allclose( ref_fv.values, fs.data_matrix[0], rtol=1e-05 )
        finally:
            rmtree( tempdir )

        from os import mkdir
        toptempdir = mkdtemp()
        try:
            class_names = []
            for letter in 'CBA':
                dirname = toptempdir + sep + letter
                mkdir( dirname )
                copy( orig_img_filepath, dirname )

            fs = FeatureSpace.NewFromDirectory( toptempdir, quiet=False, )
            self.assertEqual( fs.class_names, ['A', 'B', 'C' ] )
            for row_of_features in fs.data_matrix:
                self.assertTrue( compare( row_of_features, ref_fv.values ) )

        finally:
            rmtree( toptempdir )

#    def test_Load_GroundTruthLabels_and_Values( self ):
#        """For continuous data, we expect a float ground truth value for every sample.
#        For categorized (classed, discrete) data we expect a string label with optional
#        ground truth float for each sample. Check to make sure the ground truth makes it into
#        the right place.
#        
#        test-l.fit has class names 2cell and 4cell, which means wndcharm should try to
#        pull the 2 and the 4 from the class names and have those be the values."""
#        pass

    # --------------------------------------------------------------------------
    #@unittest.skip('')
    def test_Normalize( self ):
	"""Load unnormalized feature space, normalize,
        then compare to stored normalized feature space."""

        raw_fs = FeatureSpace.NewFromFitFile( self.test_fit_path ) 
        result_fs = raw_fs.Normalize( inplace=False )
        target_fs = FeatureSpace.NewFromFitFile( self.test_normalized_fit_path )

        from numpy.testing import assert_allclose
        assert_allclose( result_fs.data_matrix, target_fs.data_matrix, rtol=1e-05 )

        # Create reference Z-score feature space
        from scipy.stats.mstats import zscore
        from wndcharm.utils import ReplaceNonReal
        ReplaceNonReal( raw_fs.data_matrix )
        oldsettings = np.seterr( all='ignore' )
        target_fs = zscore( raw_fs.data_matrix )
        target_fs[ np.isnan( target_fs) ] = 0
        np.seterr( **oldsettings )

        result_fs = raw_fs.Normalize( inplace=False, zscore=True )
        assert_allclose( result_fs.data_matrix, target_fs )

        # See the problem with using all close...?
        # AssertionError: 
        # Not equal to tolerance rtol=1e-07, atol=0
        # (mismatch 100.0%)
        # x: array([[  18.540434,   44.441657,   30.894861, ...,    0.      ,
        #          56.162296,   48.702817],
        #       [   0.      ,    0.      ,   63.576585, ...,    0.669004,...
        # y: array([[  18.5404  ,   44.4417  ,   30.8949  , ...,    0.      ,
        #          56.1623  ,   48.7028  ],
        #       [   0.      ,    0.      ,   63.5766  , ...,    0.669004,...


        #from wndcharm.utils import compare
        #if not compare( result_fs.data_matrix, target_fs.data_matrix ):
        #    # actually, the one good thing about "allclose" is the error reporting
        #    from numpy.testing import assert_allclose
        #    assert_allclose( result_fs.data_matrix, target_fs.data_matrix )


    # --------------------------------------------------------------------------
    #@unittest.skip('')
    def test_LDATransform( self ):
	"""LDA transform"""

        tempdir = mkdtemp()
        import zipfile
        zf = zipfile.ZipFile( pychrm_test_dir + sep + 'lymphoma_iicbu2008_subset_EOSIN_ONLY_t5x6_v3.2features.fit.zip', mode='r')
        zf.extractall( tempdir )
        fitfile_path = tempdir + sep + 'lymphoma_iicbu2008_subset_eosin_t5x6_v3.2features.fit'
        try:
            kwargs = {}
            kwargs['pathname'] = fitfile_path
            kwargs['quiet'] = True
            # sampling opts: -l -t5x6 implies 5 columns and 6 rows ... I know it's weird.
            kwargs['long'] = True
            kwargs['tile_num_rows'] = 6
            kwargs['tile_num_cols'] = 5

            # against self:
            fs = FeatureSpace.NewFromFitFile( **kwargs )
            self_transformed = fs.LDATransform( reference_features=None, inplace=False )

            fit_on_fit_LDA_result = FeatureSpaceClassification.NewWND5(
                    self_transformed, self_transformed, feature_weights=None )

            # against other:
            train, test = fs.Split()
            train.LDATransform( reference_features=None, inplace=True )
            test.LDATransform( reference_features=train, inplace=True )

            split_LDA_result = FeatureSpaceClassification.NewWND5(
                    train, test, feature_weights=None )
        finally:
            rmtree( tempdir )

    # --------------------------------------------------------------------------
    #@unittest.skip('')
    def test_ClassSortSamplesByGroundTruth( self ):
        """If class names have interpolatable value, sort by that value,
        otherwise by class name alphabetical order"""

        fs = FeatureSpace( name='training set', num_samples=4, num_features=5)
        for i in (1,2,3,4,):
            fs.data_matrix[i-1].fill(i)
        fs.feature_names = [ "feat" + char for char in "ABCDE" ]
        fs._contiguous_sample_group_ids = range(4)
        fs._contiguous_sample_sequence_ids = [1]*4
        fs._contiguous_sample_names = ['ones', 'twos', 'threes', 'fours']
        fs._contiguous_ground_truth_labels = ['ones_1.0', 'twos_2.0', 'threes_3.0', 'fours_4.0']
        fs._contiguous_ground_truth_values = [1,2,3,4,]
        fs.num_classes = 4
        fs.class_names = ['ones_1.0', 'twos_2.0', 'threes_3.0', 'fours_4.0']
        fs.interpolation_coefficients = [1.0,2.0,3.0,4.0]
        fs.SortSamplesByGroundTruth(inplace=True)

        # i.e., not ['fours', 'ones', 'threes', 'twos']
        # aka, class alphabetical order
        self.assertEqual( fs.class_names, ['ones_1.0', 'twos_2.0', 'threes_3.0', 'fours_4.0'] )
        fs.SortSamplesByGroundTruth(inplace=True, force_use_labels=True)
        self.assertEqual( fs.class_names, ['fours_4.0', 'ones_1.0', 'threes_3.0', 'twos_2.0'] )



if __name__ == '__main__':
    unittest.main()
