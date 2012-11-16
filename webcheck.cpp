

//////////////////////////////////////////////////////////////////
// Borrowed from Plink v0.99s                                   //
//                                                              //
//////////////////////////////////////////////////////////////////


#define  PORT_NUM                80     
#define  IP_ADDR    "160.129.37.40"  
#define  GET_STRING "GET /wasp/files/version.txt HTTP/1.1\nHost: chgr.mc.vanderbilt.edu\nConnection: close\n\n"


void webcheck(vector<string> a)
{

#ifdef SKIP
	opts::printLog("Web-check not implemented on this system...\n");
  return;
#else
  
  opts::printLog("Web-based version check ( --noweb to skip )\n");
  
  vector<string> tokens = socketConnection( 
					    IP_ADDR,
					    PORT_NUM,
					    GET_STRING);
					    
  bool print = false;
  bool print2 = false;
  bool version_okay = true;
  vector<string>::iterator fiter;

  for (int i=0; i<tokens.size(); i++)
    {

      if (tokens[i]=="END") break;

      if (tokens[i]=="END-MESSAGE")
	{
	  print2=false;
	  continue;
	}

      if (tokens[i]=="WARN")
	{
	  if ( i < tokens.size()-1 ) 
	    {
	      i++;
		  fiter = find(a.begin(), a.end(), tokens[i]);
	      if ( fiter != a.end())
		{
		  opts::printLog("\n*** ALERT ***\n*** A warning flag has been set for: "+tokens[i]+
			   "\n*** See http://chgr.mc.vanderbilt.edu/wasp/\n");
		  warnings = true;
		}
	    }
	  continue;
	}
      
      
      if (tokens[i]=="FATAL")
	{
	  if ( i < tokens.size()-1 ) 
	    {
	      i++;
		  fiter = find(a.begin(), a.end(), tokens[i]);
	      if ( fiter != a.end())
		error("A serious warning flag has been set for: "+tokens[i]+
		    "\nWasp has been instructed to stop"+
 	            "\nPlease see http://chgr.mc.vanderbilt.edu/wasp/\n");
	    }
	  continue;
	}


      if (tokens[i]=="MESSAGE-ALL")
	{
	  print2=true;
	  continue;
	}

      // Display any other messages
      // Either conditional on old version (print)
      // or a broadcast to all users (print2)

      if ( ( print && !version_okay) || print2 ) 
	{
	  if (tokens[i]=="\\n")
	    opts::printLog("\n");
	  else
	    opts::printLog(tokens[i]+" ");
	}

      // Check version code
      if (tokens[i]=="WASPVER") 
	{
	  print=true;
	  if ( i < tokens.size() - 1) 
	    {
	      if (tokens[i+1] == _WASPVER_)
		opts::printLog(" OK, v"+_WASPVER_+" is current\n");
	      else
		{
		  opts::printLog("\n\n          *** UPDATE REQUIRED ***\n\n");
		  opts::printLog("\tThis version        : "+_WASPVER_+"\n");
		  opts::printLog("\tMost recent version : "+tokens[i+1]+"\n\n");
		  opts::printLog("Please upgrade your version of Wasp as soon as possible!\n"); 
		  opts::printLog("  (visit the above website for free download)\n\n");
		  version_okay=false;
		}

	      // Skip the version number
	      i++;
	    }
	}

    }

  // did we get the information we needed?
  if (!print) opts::printLog(" problem connecting to web\n");

  opts::printLog("\n");
#endif  
}

