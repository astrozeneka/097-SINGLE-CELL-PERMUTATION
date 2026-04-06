
export function useSshCredentials() {
    const privateKey = process.env.NEXT_PUBLIC_PRIVATE_KEY || '';
    const linuxUser = process.env.NEXT_PUBLIC_LINUX_USER || 'john';

    const credentials = {
        linux_user: linuxUser,
        private_key: privateKey.replace(/\\n/g, '\n')
    };

    return credentials;
}