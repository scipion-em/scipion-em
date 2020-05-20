import paramiko
import os

class Connect():
    def __init__(self, host, port,
                 username, password,
                 keyfilepath, keyfiletype,
                 remote_path, projectName):
        """
        create_sftp_client(host, port, username, password, keyfilepath, keyfiletype) -> SFTPClient

        Creates a SFTP client connected to the supplied host on the supplied port authenticating as the user with
        supplied username and supplied password or with the private key in a file with the supplied path.
        If a private key is used for authentication, the type of the keyfile needs to be specified as DSA or RSA.
        :rtype: SFTPClient object.

        remote_path: all paths are relative to this directory
        """
        self.sftp = None
        key = None
        self.transport = None
        try:
            if keyfilepath is not None:
                # Get private key used to authenticate user.
                if keyfiletype == 'DSA':
                    # The private key is a DSA type key.
                    key = paramiko.DSSKey.from_private_key_file(keyfilepath)
                else:
                    # The private key is a RSA type key.
                    key = paramiko.RSAKey.from_private_key_file(keyfilepath)

            # Create Transport object using supplied method of authentication

            self.transport = paramiko.Transport((host, port))
            self.transport.connect(None, username, password, key)

            self.sftp = paramiko.SFTPClient.from_transport(self.transport)

        except Exception as e:
            print('An error occurred creating SFTP client: %s: %s' % (e.__class__, e))
            if self.sftp is not None:
                self.sftp.close()
            if self.transport is not None:
                self.transport.close()
        self.remote_path = remote_path
        directory = os.path.join(self.remote_path, projectName)
        try:
            self.sftp.chdir(directory)  # Test if directory exists
        except IOError:
            self.sftp.mkdir(directory)  # Create directory

    def put(self,listLocalPaths, listRemotePaths):
        try:
            for local, remote in zip(listLocalPaths, listRemotePaths):
                remote = os.path.join(self.remote_path, remote)
                print(local, "-->", remote)
                self.sftp.put(local, remote, confirm=True)
        except IOError:
            pass
        except OSError:
            pass
        except Exception as e:
            print(str(e))
        return 0

    def close(self):
        if self.sftp is not None:
            self.sftp.close()
        if self.transport is not None:
            self.transport.close()


